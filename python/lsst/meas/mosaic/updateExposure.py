#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import re
import numpy

from . import getFCorImg, FluxFitParams, getJImg, calculateJacobian
import lsst.geom as geom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.pipe.base import Struct, TaskError
from lsst.daf.persistence import NoResults
from . import utils as mosaicUtils

__all__ = ("applyMosaicResults", "getMosaicResults", "applyMosaicResultsExposure",
           "applyMosaicResultsCatalog", "applyCalib")


def applyMosaicResults(dataRef, calexp=None):
    """Deprecated function to apply the results to an exposure

    Deprecated, because the mosaic results can be applied to more than
    one kind of target, so it's worth changing the name to be specific.
    """
    return applyMosaicResultsExposure(dataRef, calexp).exposure


def applyMosaicResultsExposure(dataRef, calexp=None):
    """Update an Exposure with the Wcs, Calib, and flux scaling from meas_mosaic.

    If None, the calexp will be loaded from the dataRef.  Otherwise it is
    updated in-place.

    This assumes that the mosaic solution exists; an exception will be raised
    in the event that it does not.
    """
    if calexp is None:
        calexp = dataRef.get("calexp", immediate=True)

    nQuarter = calexp.getDetector().getOrientation().getNQuarter()
    dims = calexp.getDimensions()
    hscRun = mosaicUtils.checkHscStack(calexp.getMetadata())

    # Need the dimensions in coordinates used by meas_mosaic which defines 0,0 as the
    # lower-left hand corner on the sky
    if hscRun is None:
        if nQuarter%2 != 0:
            width, height = calexp.getDimensions()
            dims = geom.Extent2I(height, width)

    # return results in meas_mosaic coordinate system
    mosaic = getMosaicResults(dataRef, dims)

    # rotate wcs back to LSST coordinate system
    if nQuarter%4 != 0 and hscRun is None:
        import lsst.meas.astrom as measAstrom
        mosaic.wcs = measAstrom.rotateWcsPixelsBy90(mosaic.wcs, 4 - nQuarter, dims)
    calexp.setWcs(mosaic.wcs)

    fluxMag0 = mosaic.calib.getInstFluxAtZeroMagnitude()
    calexp.setPhotoCalib(afwImage.makePhotoCalibFromCalibZeroPoint(fluxMag0, 0.0))

    mi = calexp.getMaskedImage()
    # rotate photometric correction to LSST coordiantes
    if nQuarter%4 != 0 and hscRun is None:
        mosaic.fcor = afwMath.rotateImageBy90(mosaic.fcor, 4 - nQuarter)
    mi *= mosaic.fcor

    return Struct(exposure=calexp, mosaic=mosaic)


def getFluxFitParams(dataRef):
    """Retrieve the flux correction parameters determined by meas_mosaic

    If the flux correction parameters do not exist, an exception will
    be raised.
    """
    calexp_md = dataRef.get("calexp_md", immediate=True)
    hscRun = mosaicUtils.checkHscStack(calexp_md)
    if hscRun is not None:
        ffpHeader = dataRef.get("fcr_hsc_md", immediate=True)
    else:
        ffpHeader = dataRef.get("fcr_md", immediate=True)
    photoCalib = dataRef.get("fcr_photoCalib")
    ffp = FluxFitParams(ffpHeader)

    wcs = getWcs(dataRef)

    if hscRun is None:
        detector = dataRef.get("camera")[dataRef.dataId["ccd"]]
        nQuarter = detector.getOrientation().getNQuarter()
        if nQuarter%4 != 0:
            # Have to put this import here due to circular dependence in forcedPhotCcd.py in meas_base
            import lsst.meas.astrom as measAstrom
            dimensions = dataRef.get("calexp_bbox").getDimensions()
            wcs = measAstrom.rotateWcsPixelsBy90(wcs, nQuarter, dimensions)
    return Struct(ffp=ffp, calib=photoCalib, wcs=wcs)


def getWcs(dataRef):
    """Retrieve the Wcs determined by meas_mosaic

    If the Wcs does not exist, an exception will be raised.
    """
    calexp_md = dataRef.get("calexp_md", immediate=True)
    hscRun = mosaicUtils.checkHscStack(calexp_md)
    if hscRun is not None:
        # Backwards compatibility with the very oldest meas_mosaic outputs
        return dataRef.get("wcs_hsc").getWcs()
    try:
        # Modern meas_mosaic outputs.
        return dataRef.get("jointcal_wcs")
    except NoResults:
        # Backwards compatibility with old meas_mosaic outputs
        return dataRef.get("wcs").getWcs()


def getMosaicResults(dataRef, dims=None):
    """Retrieve the results of meas_mosaic

    If None, the dims will be determined from the calexp header.
    """
    ffp = getFluxFitParams(dataRef)

    if dims is None:
        bbox = dataRef.get("calexp_bbox")
        width, height = bbox.getWidth(), bbox.getHeight()
    else:
        width, height = dims

    fcor = getFCorImg(ffp.ffp, width, height)
    jcor = getJImg(ffp.wcs, width, height)
    fcor *= jcor
    del jcor

    return Struct(wcs=ffp.wcs, calib=ffp.calib, fcor=fcor)


def applyMosaicResultsCatalog(dataRef, catalog, addCorrection=True):
    """!Apply the results of meas_mosaic to a source catalog

    The coordinates and all fluxes are updated in-place with the meas_mosaic solution.

    This assumes that the mosaic solution exists; an exception will be raised
    in the event that it does not.
    """
    ffp = getFluxFitParams(dataRef)
    calexp_md = dataRef.get("calexp_md", immediate=True)
    hscRun = mosaicUtils.checkHscStack(calexp_md)
    if hscRun is None:
        detector = dataRef.get("camera")[dataRef.dataId["ccd"]]
        nQuarter = detector.getOrientation().getNQuarter()
        if nQuarter%4 != 0:
            dimensions = dataRef.get("calexp_bbox").getDimensions()
            catalog = mosaicUtils.rotatePixelCoords(catalog, dimensions.getX(), dimensions.getY(),
                                                    nQuarter)
    xx, yy = catalog.getX(), catalog.getY()
    corr = numpy.power(10.0, -0.4*ffp.ffp.eval(xx, yy))*calculateJacobian(ffp.wcs, xx, yy)

    if addCorrection:
        mapper = afwTable.SchemaMapper(catalog.schema, True)
        for s in catalog.schema:
            mapper.addMapping(s.key)
        corrField = afwTable.Field[float]("mosaic_corr", "Magnitude correction from meas_mosaic")
        corrKey = mapper.addOutputField(corrField)
        outCatalog = type(catalog)(mapper.getOutputSchema())
        outCatalog.extend(catalog, mapper=mapper)
        outCatalog[corrKey][:] = corr
        catalog = outCatalog

    fluxKeys, errKeys = getFluxKeys(catalog.schema, hscRun=hscRun)
    for name, key in list(fluxKeys.items()) + list(errKeys.items()):
        # Note this skips correcting the aperture fluxes in HSC processed data, but that's ok because
        # we are using the flux_sinc as our comparison to base_CircularApertureFlux_12_0_flux
        if key.subfields is None:
            catalog[key][:] *= corr

    # Now rotate them back to the LSST coord system
    if hscRun is None:
        if nQuarter%4 != 0:
            catalog = mosaicUtils.rotatePixelCoordsBack(catalog, dimensions.getX(),
                                                        dimensions.getY(), nQuarter)

    wcs = getWcs(dataRef)
    for rec in catalog:
        rec.updateCoord(wcs)

    return Struct(catalog=catalog, wcs=wcs, ffp=ffp)


def applyCalib(catalog, photoCalib, hscRun=None):
    """Convert all fluxes in a catalog to magnitudes

    The fluxes are converted in-place, so that the "_flux*" are now really
    magnitudes.
    """
    fluxKeys, errKeys = getFluxKeys(catalog.schema, hscRun=hscRun)
    mapper = afwTable.SchemaMapper(catalog.schema, True)
    for item in catalog.schema:
        name = item.field.getName()
        if name in fluxKeys:
            continue
        mapper.addMapping(item.key)
    aliasMap = catalog.schema.getAliasMap()

    newFluxKeys = {}
    newErrKeys = {}
    for name in fluxKeys:
        fluxField = catalog.schema.find(name).field
        newName = name.replace("instFlux", "mag")
        newField = fluxField.__class__(newName, "Calibrated magnitude from %s (%s)" %
                                       (fluxField.getName(), fluxField.getDoc()), "mag")
        newFluxKeys[newName] = mapper.addMapping(fluxKeys[name], newField)

        errName = "Err"
        if hscRun is not None:
            errName = "_err"

        if name + errName in errKeys:
            errField = catalog.schema.find(name + errName).field
            newErrField = errField.__class__(newName + errName,
                                             "Calibrated magnitude error from %s (%s)" %
                                             (errField.getName(), errField.getDoc()), "mag")
            newErrKeys[newName] = mapper.addMapping(errKeys[name + errName], newErrField)
        aliasMap.set(name, newName)
        aliasMap.set(name + errName, newName + errName)

    newCatalog = afwTable.SourceCatalog(mapper.getOutputSchema())
    newCatalog.extend(catalog, mapper=mapper)

    for name, key in newFluxKeys.items():
        flux = newCatalog[key]
        if name in newErrKeys:
            result = photoCalib.instFluxToMagnitude(newCatalog, name.strip('_mag'))
            flux[:] = result[:, 0]
            newCatalog[newErrKeys[name]] = result[:, 1]
        else:
            flux[:] = numpy.array([photoCalib.instFluxToMagnitude(f) for f in flux])

    return newCatalog


def getFluxKeys(schema, hscRun=None):
    """Retrieve the flux and flux error keys from a schema

    Both are returned as dicts indexed on the flux name (e.g. "base_PsfFlux" or "base_CmodelFlux").
    """
    if hscRun is None:
        fluxTypeStr = "_instFlux"
        fluxSchemaItems = schema.extract("*" + fluxTypeStr)
        # Do not include any flag fields (as determined by their type).  Also exclude
        # slot fields, as these would effectively duplicate whatever they point to.
        fluxKeys = dict((name, schemaItem.key) for name, schemaItem in list(fluxSchemaItems.items()) if
                        schemaItem.field.getTypeString() != "Flag" and
                        not name.startswith("slot"))

        errSchemaItems = schema.extract("*" + fluxTypeStr + "Err")
        errKeys = dict((name, schemaItem.key) for name, schemaItem in list(errSchemaItems.items()) if
                       name[:-len("Err")] in fluxKeys)
    else:
        schemaKeys = dict((s.field.getName(), s.key) for s in schema)
        fluxKeys = dict((name, key) for name, key in schemaKeys.items() if
                        re.search(r"^(flux\_\w+|\w+\_instFlux)$", name) and not
                        re.search(r"^(\w+\_apcorr)$", name) and name + "_err" in schemaKeys)
        errKeys = dict((name + "_err", schemaKeys[name + "_err"]) for name in fluxKeys if
                       name + "_err" in schemaKeys)

    if len(fluxKeys) == 0:
        raise TaskError("No flux keys found")

    return fluxKeys, errKeys
