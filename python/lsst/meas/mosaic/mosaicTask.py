#!/usr/bin/env python

import os
import math
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.cameraGeom.utils        as cameraGeomUtils
import lsst.afw.geom                    as afwGeom
import lsst.afw.image                   as afwImage
import lsst.afw.table                   as afwTable
import lsst.pex.config                  as pexConfig
import lsst.pipe.base                   as pipeBase
import lsst.meas.mosaic.mosaicLib       as measMosaic
import lsst.meas.astrom.astrom          as measAstrom
from lsst.meas.photocal.colorterms import Colorterm
import hsc.pipe.base.camera             as hscCamera

class MosaicConfig(pexConfig.Config):
    nBrightest = pexConfig.RangeField(
        doc="number of stars used for fitting per exposure",
        dtype=int,
        default=300, min=100)
    radXMatch = pexConfig.RangeField(
        doc="radius to cross-match objects between expsoures in arcsec",
        dtype=float,
        default=5.0, min=3.0)
    fittingOrder = pexConfig.RangeField(
        doc="fitting order",
        dtype=int,
        default=5, min=2)
    internalFitting = pexConfig.Field(
        doc="Use stars without catalog matching for fitting?",
        dtype=bool,
        default=True)
    solveCcd = pexConfig.Field(
        doc="Solve CCD alignment?",
        dtype=bool,
        default=True)
    allowRotation = pexConfig.Field(
        doc="Solve rotation?",
        dtype=bool,
        default=True)
    catRMS = pexConfig.Field(
        doc="Positional error in reference catalog (degree)",
        dtype=float,
        default=0.040/3600.)
    chebyshev = pexConfig.Field(
        doc="Use Chebyshev polynomials for flux fitting?",
        dtype=bool,
        default=True)
    fluxFitOrder = pexConfig.RangeField(
        doc="flux fitting order",
        dtype=int,
        default=5, min=0)
    fluxFitAbsolute = pexConfig.Field(
        doc="Fit to catalog flux?",
        dtype=bool,
        default=False)
    outputDir = pexConfig.Field(
        doc="Output directory to write diagnostics plots",
        dtype=str,
        default=".")
    outputDiag = pexConfig.Field(
        doc="Output diagnostics plots",
        dtype=bool,
        default=False)

class MosaicTask(pipeBase.Task):

    ConfigClass = MosaicConfig
    _DefaultName = "meas_mosaic"

    def readCcd(self, camera, ccdIds):
        self.log.info("Reading CCD info ...")

        ccds = measMosaic.CcdSet()
        for i in ccdIds:
            ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(i)))
            ccds[i] = ccd

        # Calculate mean position of all CCD chips
        sx = sy = 0.
        for ccd in ccds.values():
            center = ccd.getCenter().getPixels(ccd.getPixelSize())
            sx += center[0]
            sy += center[1]
        dx = sx / ccds.size()
        dy = sy / ccds.size()

        # Shift the origin of CCD chips
        for ccd in ccds.values():
            pixelSize = ccd.getPixelSize()
            ccd.setCenter(ccd.getCenter()
                          - cameraGeom.FpPoint(dx*pixelSize, dy*pixelSize)
                          - cameraGeom.FpPoint(ccd.getAllPixels(True).getWidth()*0.5,
                                               ccd.getAllPixels(True).getHeight()*0.5))
        
        return ccds
        
    def getWcsForCcd(self, butler, frameId, ccdId):

        dataId = {'visit': frameId, 'ccd':ccdId}

        try:
            md = butler.get('calexp_md', dataId)
            return afwImage.makeWcs(md)
        except Exception, e:
            print "Failed to read: %s for %s" % (e, dataId)
            return None

    def readWcs(self, butler, frameIds, ccdSet):
        
        self.log.info("Reading WCS ...")

        wcsDic = measMosaic.WcsDic()
        for frameId in frameIds:
            found = 0
            for ccdId in ccdSet.keys():
                ccd = ccdSet[ccdId]
                dataId = {'visit': frameId, 'ccd': ccdId}
                if (butler.datasetExists('calexp',  dataId) and
                    butler.datasetExists('src',     dataId) and
                    butler.datasetExists('icSrc',   dataId) and
                    butler.datasetExists('icMatch', dataId)):
                    found = 1
                    break
            if found:
                wcs = self.getWcsForCcd(butler, frameId, ccdId)
                ccd = ccdSet[ccdId]
                offset = ccd.getCenter().getPixels(ccd.getPixelSize())
                wcs.shiftReferencePixel(offset[0], offset[1])
                wcsDic[frameId] = wcs

        return wcsDic

    def removeNonExistCcd(self, butler, ccdSet, wcsDic):
        num = dict()
        for ichip in ccdSet.keys():
            num[ichip] = 0

        for iexp in wcsDic.keys():
            for ichip in ccdSet.keys():
                dataId = {'visit': iexp, 'ccd': ichip}
                if (butler.datasetExists('calexp',  dataId) and
                    butler.datasetExists('src',     dataId) and
                    butler.datasetExists('icSrc',   dataId) and
                    butler.datasetExists('icMatch', dataId)):
                    num[ichip] += 1

        for ichip in ccdSet.keys():
            if num[ichip] == 0:
                ccdSet.erase(ichip)
            
    def selectStars(self, sources, includeSaturated=False):
        if len(sources) == 0:
            return []
        if isinstance(sources, afwTable.SourceCatalog):
            extended = sources.columns["classification.extendedness"]
            saturated = sources.columns["flags.pixel.saturated.center"]
            try:
                nchild = sources.columns["deblend.nchild"]
            except:
                nchild = numpy.zeros(len(sources))
            indices = numpy.where(numpy.logical_and(numpy.logical_and(extended < 0.5, saturated == False), nchild == 0))[0]
            return [sources[int(i)] for i in indices]

        psfKey = None                       # Table key for classification.psfstar
        if isinstance(sources, afwTable.ReferenceMatchVector) or isinstance(sources[0], afwTable.ReferenceMatch):
            sourceList = [s.second for s in sources]
            psfKey = sourceList[0].schema.find("classification.psfstar").getKey()
        else:
            sourceList = sources

        schema = sourceList[0].schema
        extKey = schema.find("classification.extendedness").getKey()
        satKey = schema.find("flags.pixel.saturated.center").getKey()

        stars = []
        for includeSource, checkSource in zip(sources, sourceList):
            star = (psfKey is not None and checkSource.get(psfKey)) or checkSource.get(extKey) < 0.5
            saturated = checkSource.get(satKey)
            if star and (includeSaturated or not saturated):
                stars.append(includeSource)
        return stars

    def getAllForCcd(self, butler, astrom, frame, ccd, ct=None):

        data = {'visit': frame, 'ccd': ccd}

        try:
            if not butler.datasetExists('src', data):
                raise RuntimeError("no data for src %s" % (data))
            if not butler.datasetExists('calexp_md', data):
                raise RuntimeError("no data for calexp_md %s" % (data))
            md = butler.get('calexp_md', data)
            wcs = afwImage.makeWcs(md)

            sources = butler.get('src', data)
            if False:
                matches = measAstrom.readMatches(butler, data)
            else:
                icSrces = butler.get('icSrc', data)
                packedMatches = butler.get('icMatch', data)
                matches = astrom.joinMatchListWithCatalog(packedMatches, icSrces, True)
                if ct != None:
                    if matches[0].first != None:
                        refSchema = matches[0].first.schema
                    else:
                        refSchema = matches[1].first.schema
                    key_p = refSchema.find(ct.primary).key
                    key_s = refSchema.find(ct.secondary).key
                    key_f = refSchema.find("flux").key
                    for m in matches:
                        if m.first != None:
                            refFlux1 = m.first.get(key_p)
                            refFlux2 = m.first.get(key_s)
                            refMag1 = -2.5*math.log10(refFlux1)
                            refMag2 = -2.5*math.log10(refFlux2)
                            refMag = ct.transformMags(ct.primary, refMag1, refMag2)
                            refFlux = math.pow(10.0, -0.4*refMag)
                            if refFlux == refFlux:
                                m.first.set(key_f, refFlux)
                            else:
                                m.first = None

            sources = self.selectStars(sources)
            selMatches = self.selectStars(matches)
            if len(selMatches) < 10:
                matches = self.selectStars(matches, True)
            else:
                matches = selMatches
        except Exception, e:
            print "Failed to read: %s" % (e)
            return None, None, None
    
        return sources, matches, wcs

    def readCatalog(self, butler, frameIds, ccdIds, ct=None):
        self.log.info("Reading catalogs ...")

        sourceSet = measMosaic.SourceGroup()
        matchList = measMosaic.SourceMatchGroup()
        astrom = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
        for frameId in frameIds:
            ss = []
            ml = []
            for ccdId in ccdIds:
                sources, matches, wcs = self.getAllForCcd(butler, astrom, frameId, ccdId, ct)
                if sources != None:
                    for s in sources:
                        if numpy.isfinite(s.getRa().asDegrees()): # get rid of NaN
                            src = measMosaic.Source(s)
                            src.setExp(frameId)
                            src.setChip(ccdId)
                            ss.append(src)
                    for m in matches:
                        if m.first != None and m.second != None:
                            match = measMosaic.SourceMatch(measMosaic.Source(m.first, wcs), measMosaic.Source(m.second))
                            match.second.setExp(frameId)
                            match.second.setChip(ccdId)
                            ml.append(match)
            sourceSet.push_back(ss)
            matchList.push_back(ml)

        return sourceSet, matchList

    def countObsInSourceGroup(self, sg):
        num = 0
        for s in sg:
            num += (len(s) - 1)

        return num

    def mergeCatalog(self, sourceSet, matchList, ccdSet, d_lim, nbrightest):

        self.log.info("Creating kd-tree for matched catalog ...")
        self.log.info("len(matchList) = "+str(len(matchList))+" "+
                      str([len(matches) for matches in matchList]))
        rootMat = measMosaic.kdtreeMat(matchList)
        allMat = rootMat.mergeMat()
        self.log.info("# of allMat : %d" % self.countObsInSourceGroup(allMat))
        self.log.info('len(allMat) = %d' % len(allMat))
    
        self.log.info("Creating kd-tree for source catalog ...")
        self.log.info('len(sourceSet) = '+str(len(sourceSet))+" "+
                      str([len(sources) for sources in sourceSet]))
        rootSource = measMosaic.kdtreeSource(sourceSet, rootMat, ccdSet, d_lim, nbrightest)
        allSource = rootSource.mergeSource()
        self.log.info("# of allSource : %d" % self.countObsInSourceGroup(allSource))
        self.log.info('len(allSource) = %d' % len(allSource))

        return allMat, allSource

    def writeNewWcs(self):
        self.log.info("Write New WCS ...")
        exp = afwImage.ExposureI(0,0)
        for iexp in self.coeffSet.keys():
            for ichip in self.ccdSet.keys():
                c = measMosaic.convertCoeff(self.coeffSet[iexp], self.ccdSet[ichip]);
                wcs = measMosaic.wcsFromCoeff(c);
                exp.setWcs(wcs)
                scale = self.fexp[iexp] * self.fchip[ichip]
                calib = afwImage.Calib()
                calib.setFluxMag0(1.0/scale)
                exp.setCalib(calib)
                try:
                    self.butler.put(exp, 'wcs', dict(visit=iexp, ccd=ichip))
                except Exception, e:
                    print "failed to write something: %s" % (e)

    def writeFcr(self):
        self.log.info("Write Fcr ...")
        for iexp in self.coeffSet.keys():
            for ichip in self.ccdSet.keys():
                newP = measMosaic.convertFluxFitParams(self.coeffSet[iexp], self.ccdSet[ichip],
                                                       measMosaic.FluxFitParams(self.ffp))
                metadata = measMosaic.metadataFromFluxFitParams(newP)
                exp = afwImage.ExposureI(0,0)
                exp.getMetadata().combine(metadata)
                scale = self.fexp[iexp] * self.fchip[ichip]
                calib = afwImage.Calib()
                calib.setFluxMag0(1.0/scale)
                exp.setCalib(calib)
                try:
                    self.butler.put(exp, 'fcr', dict(visit=iexp, ccd=ichip))
                except Exception, e:
                    print "failed to write something: %s" % (e)

    def getExtent(self, matchVec):
        u_max = float("-inf")
        v_max = float("-inf")
        for m in matchVec:
            if (math.fabs(m.u) > u_max):
                u_max = math.fabs(m.u)
            if (math.fabs(m.v) > v_max):
                v_max = math.fabs(m.v)

        return u_max, v_max

    def checkInputs(self, wcsDic, sourceSet, matchList):
        newWcsDic = measMosaic.WcsDic()
        newSourceSet = measMosaic.SourceGroup()
        newMatchList = measMosaic.SourceMatchGroup()
        for i, (wcs, frame, ss, ml) in enumerate(zip(wcsDic.values(), wcsDic.keys(), sourceSet, matchList)):
            if len(ss) > 0 or len(ml) > 0:
                newWcsDic[frame] = wcs
                newSourceSet.push_back(ss)
                newMatchList.push_back(ml)
        return newWcsDic, newSourceSet, newMatchList

    def plotCcd(self, coeffx0, coeffy0):
        for ccd in self.ccdSet.values():
            w = ccd.getAllPixels(True).getWidth()
            h = ccd.getAllPixels(True).getHeight()
            x0 = ccd.getCenter().getPixels(ccd.getPixelSize())[0] + coeffx0
            y0 = ccd.getCenter().getPixels(ccd.getPixelSize())[1] + coeffy0
            t0 = ccd.getOrientation().getYaw()
            x = numpy.array([x0,
                             x0 + w * math.cos(t0),
                             x0 + w * math.cos(t0) - h * math.sin(t0),
                             x0 - h * math.sin(t0),
                             x0])
            y = numpy.array([y0,
                             y0 + w * math.sin(t0),
                             y0 + w * math.sin(t0) + h * math.cos(t0),
                             y0 + h * math.cos(t0),
                             y0])
            plt.plot(x, y, 'k-')

    def plotJCont(self, iexp):
        coeff = self.coeffSet[iexp]

        scale = coeff.pixelScale()
        deg2pix = 1. / scale

        delta = 300.
        if (self.ccdSet.size() >= 100):
            x = numpy.arange(-18000., 18000., delta)
            y = numpy.arange(-18000., 18000., delta)
            levels = numpy.linspace(0.81, 1.02, 36)
        else:
            x = numpy.arange(-6000., 6000., delta)
            y = numpy.arange(-6000., 6000., delta)
            levels = numpy.linspace(0.88, 1.02, 36)
        X, Y = numpy.meshgrid(x, y)
        Z = numpy.zeros((len(X),len(Y)))

        for j in range(len(Y)):
            for i in range(len(X)):
                Z[i][j] = coeff.detJ(X[i][j], Y[i][j]) * deg2pix ** 2

        plt.clf()
        plt.contourf(X, Y, Z, levels=levels)
        plt.colorbar()

        self.plotCcd(coeff.x0, coeff.y0)

        plt.savefig(os.path.join(self.outputDir, "jcont_%d.png" % (iexp)), format='png')

    def plotFCorCont(self, iexp):
        coeff = self.coeffSet[iexp]

        delta = 300.
        if (self.ccdSet.size() >= 100):
            x = numpy.arange(-18000., 18000., delta)
            y = numpy.arange(-18000., 18000., delta)
            levels = numpy.linspace(0.81, 1.02, 36)
        else:
            x = numpy.arange(-6000., 6000., delta)
            y = numpy.arange(-6000., 6000., delta)
            levels = numpy.linspace(0.86, 1.14, 36)
        X, Y = numpy.meshgrid(x, y)
        Z = numpy.zeros((len(X),len(Y)))

        for j in range(len(Y)):
            for i in range(len(X)):
                Z[i][j] = 10**(-0.4*self.ffp.eval(X[i][j], Y[i][j]))

        plt.clf()
        plt.contourf(X, Y, Z, levels=levels)
        plt.colorbar()

        self.plotCcd(coeff.x0, coeff.y0)
        
        plt.savefig(os.path.join(self.outputDir, "fcont_%d.png" % (iexp)), format='png')

    def plotResPosArrow2D(self, iexp):
        _xm = []
        _ym = []
        _dxm = []
        _dym = []
        for m in self.matchVec:
            if (m.good == True and m.iexp == iexp):
                _xm.append(m.u)
                _ym.append(m.v)
                _dxm.append((m.xi_fit - m.xi) * 3600)
                _dym.append((m.eta_fit - m.eta) * 3600)
        _xs = []
        _ys = []
        _dxs = []
        _dys = []
        if (self.sourceVec != None):
            for s in self.sourceVec:
                if (s.good == True and s.iexp == iexp):
                    _xs.append(s.u)
                    _ys.append(s.v)
                    _dxs.append((s.xi_fit - s.xi) * 3600)
                    _dys.append((s.eta_fit - s.eta) * 3600)

        xm = numpy.array(_xm)
        ym = numpy.array(_ym)
        dxm = numpy.array(_dxm)
        dym = numpy.array(_dym)
        xs = numpy.array(_xs)
        ys = numpy.array(_ys)
        dxs = numpy.array(_dxs)
        dys = numpy.array(_dys)

        plt.clf()
        plt.rc('text', usetex=True)

        q = plt.quiver(xm, ym, dxm, dym, units='inches', angles='xy', scale=1, color='green')
        plt.quiverkey(q, 0, 4500, 0.1, "0.1 arcsec", coordinates='data', color='black')
        plt.quiver(xs, ys, dxs, dys, units='inches', angles='xy', scale=1, color='red')

        self.plotCcd(self.coeffSet[iexp].x0, self.coeffSet[iexp].y0)
        plt.axes().set_aspect('equal')

        plt.savefig(os.path.join(self.outputDir, "ResPosArrow2D_%d.png" % (iexp)), format='png')

    def clippedStd(self, a, n):
        avg = a.mean()
        std = a.std()
        for i in range(n):
            b = a[numpy.fabs(a-avg) < 2.1*std]
            avg = b.mean()
            std = b.std()

        b = a[numpy.fabs(a-avg) < 2.1*std]
        avg = b.mean()
        std = b.std()
            
        return [std, avg, len(b)]

    def plotResPosScatter(self):
        _x = []
        _y = []
        _xbad = []
        _ybad = []
        _xm = []
        _ym = []
        f = open(os.path.join(self.outputDir, "dpos.dat"), "wt")
        for m in self.matchVec:
            if (m.good == True):
                _x.append((m.xi_fit - m.xi) * 3600)
                _y.append((m.eta_fit - m.eta) * 3600)
                _xm.append((m.xi_fit - m.xi) * 3600)
                _ym.append((m.eta_fit - m.eta) * 3600)
                f.write("m %f %f %f %f %f %f 1\n" % (m.xi_fit, m.eta_fit,
                                                     m.xi, m.eta, m.u, m.v))
            else:
                _xbad.append((m.xi_fit - m.xi) * 3600)
                _ybad.append((m.eta_fit - m.eta) * 3600)
                f.write("m %f %f %f %f %f %f 0\n" % (m.xi_fit, m.eta_fit,
                                                     m.xi, m.eta, m.u, m.v))
        _xs = []
        _ys = []
        if (self.sourceVec != None):
            for s in self.sourceVec:
                if (s.good == True):
                    _x.append((s.xi_fit - s.xi) * 3600)
                    _y.append((s.eta_fit - s.eta) * 3600)
                    _xs.append((s.xi_fit - s.xi) * 3600)
                    _ys.append((s.eta_fit - s.eta) * 3600)
                    f.write("s %f %f %f %f %f %f 1\n" % (s.xi_fit, s.eta_fit,
                                                         s.xi, s.eta, s.u, s.v))
                else:
                    _xbad.append((s.xi_fit - s.xi) * 3600)
                    _ybad.append((s.eta_fit - s.eta) * 3600)
                    f.write("s %f %f %f %f %f %f 0\n" % (s.xi_fit, s.eta_fit,
                                                         s.xi, s.eta, s.u, s.v))
        f.close()

        d_xi = numpy.array(_x)
        d_eta = numpy.array(_y)
        d_xi_m = numpy.array(_xm)
        d_eta_m = numpy.array(_ym)
        d_xi_s = numpy.array(_xs)
        d_eta_s = numpy.array(_ys)
        d_xi_bad = numpy.array(_xbad)
        d_eta_bad = numpy.array(_ybad)

        xi_std,  xi_mean,  xi_n  = self.clippedStd(d_xi, 2)
        eta_std, eta_mean, eta_n = self.clippedStd(d_eta, 2)
        xi_std_m,  xi_mean_m,  xi_n_m  = self.clippedStd(d_xi_m, 2)
        eta_std_m, eta_mean_m, eta_n_m = self.clippedStd(d_eta_m, 2)
        xi_std_s,  xi_mean_s,  xi_n_s  = self.clippedStd(d_xi_s, 2)
        eta_std_s, eta_mean_s, eta_n_s = self.clippedStd(d_eta_s, 2)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
        plt.plot(d_xi_bad, d_eta_bad, 'k,', markeredgewidth=0)
        plt.plot(d_xi_m, d_eta_m, 'g,', markeredgewidth=0)
        plt.plot(d_xi_s, d_eta_s, 'r,', markeredgewidth=0)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)

        plt.xlabel(r'$\Delta\xi$ (arcsec)')
        plt.ylabel(r'$\Delta\eta$ (arcsec)')

        bins = numpy.arange(-0.5, 0.5, 0.01) + 0.005

        ax = plt.subplot2grid((5,6),(0,0), colspan=4)
        if self.sourceVec != None:
            plt.hist([d_xi, d_xi_m, d_xi_s], bins=bins, normed=False, histtype='step')
        else:
            plt.hist([d_xi, d_xi_m], bins=bins, normed=False, histtype='step')
        plt.text(0.75, 0.7, r"$\sigma=$%5.3f" % (xi_std), transform=ax.transAxes, color='blue')
        plt.text(0.75, 0.5, r"$\sigma=$%5.3f" % (xi_std_m), transform=ax.transAxes, color='green')
        y = mlab.normpdf(bins, xi_mean_m, xi_std_m)
        plt.plot(bins, y*xi_n_m*0.01, 'g:')
        if self.sourceVec != None:
            plt.text(0.75, 0.3, r"$\sigma=$%5.3f" % (xi_std_s), transform=ax.transAxes, color='red')
            y = mlab.normpdf(bins, xi_mean_s, xi_std_s)
            plt.plot(bins, y*xi_n_s*0.01, 'r:')
        plt.xlim(-0.5, 0.5)

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(d_eta, bins=bins, normed=False, orientation='horizontal', histtype='step')
        plt.hist(d_eta_m, bins=bins, normed=False, orientation='horizontal', histtype='step')
        if self.sourceVec != None:
            plt.hist(d_eta_s, bins=bins, normed=False, orientation='horizontal', histtype='step')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (eta_std), rotation=270, transform=ax.transAxes, color='blue')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (eta_std_m), rotation=270, transform=ax.transAxes, color='green')
        y = mlab.normpdf(bins, eta_mean_m, eta_std_m)
        plt.plot(y*eta_n_m*0.01, bins, 'g:')
        if self.sourceVec != None:
            plt.text(0.3, 0.25, r"$\sigma=$%5.3f" % (eta_std_s), rotation=270, transform=ax.transAxes, color='red')
            y = mlab.normpdf(bins, eta_mean_s, eta_std_s)
            plt.plot(y*eta_n_s*0.01, bins, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.5, 0.5)

        plt.savefig(os.path.join(self.outputDir, "ResPosScatter.png"), format='png')

    def plotMdM(self):
        _dmag_m = []
        _dmag_cat_m = []
        _dmag_s = []
        _dmag_a = []
        _dmag_bad = []
        _dmag_cat_bad = []
        _mag0_m = []
        _mag_cat_m = []
        _mag0_s = []
        _mag0_bad = []
        _mag_cat_bad = []
        f = open(os.path.join(self.outputDir, 'dmag.dat'), 'wt')
        for m in self.matchVec:
            if (m.good == True and m.mag != -9999 and m.jstar != -1 and m.mag0 != -9999 and m.mag_cat != -9999):
                mag = m.mag
                mag0 = m.mag0
                mag_cat = m.mag_cat
                exp_cor = -2.5 * math.log10(self.fexp[m.iexp])
                chip_cor = -2.5 * math.log10(self.fchip[m.ichip])
                gain_cor = self.ffp.eval(m.u, m.v)
                mag_cor = mag + exp_cor + chip_cor + gain_cor
                diff = mag_cor - mag0
                _dmag_m.append(diff)
                _dmag_a.append(diff)
                _mag0_m.append(mag0)
                _dmag_cat_m.append(mag_cor - mag_cat)
                _mag_cat_m.append(mag_cat)
                f.write("m %f %f %f %f %f 1\n" % (mag_cor, mag0, mag_cat,
                                                  m.u, m.v))
            else:
                mag = m.mag
                mag0 = m.mag0
                mag_cat = m.mag_cat
                exp_cor = -2.5 * math.log10(self.fexp[m.iexp])
                chip_cor = -2.5 * math.log10(self.fchip[m.ichip])
                gain_cor = self.ffp.eval(m.u, m.v)
                mag_cor = mag + exp_cor + chip_cor + gain_cor
                diff = mag_cor - mag0
                _dmag_bad.append(diff)
                _mag0_bad.append(mag0)
                _dmag_cat_bad.append(mag_cor - mag_cat)
                _mag_cat_bad.append(mag_cat)
                f.write("m %f %f %f %f %f 0\n" % (mag_cor, mag0, mag_cat,
                                                  m.u, m.v))
        if self.sourceVec != None:
            for s in self.sourceVec:
                if (s.good == True and s.mag != -9999 and s.jstar != -1):
                    mag = s.mag
                    mag0 = s.mag0
                    exp_cor = -2.5 * math.log10(self.fexp[s.iexp])
                    chip_cor = -2.5 * math.log10(self.fchip[s.ichip])
                    gain_cor = self.ffp.eval(s.u, s.v)
                    mag_cor = mag + exp_cor + chip_cor + gain_cor
                    diff = mag_cor - mag0
                    _dmag_s.append(diff)
                    _dmag_a.append(diff)
                    _mag0_s.append(mag0)
                    f.write("s %f %f %f %f %f 1\n" % (mag_cor, mag0, -9999,
                                                      s.u, s.v))
                else:
                    mag = s.mag
                    mag0 = s.mag0
                    exp_cor = -2.5 * math.log10(self.fexp[s.iexp])
                    chip_cor = -2.5 * math.log10(self.fchip[s.ichip])
                    gain_cor = self.ffp.eval(s.u, s.v)
                    mag_cor = mag + exp_cor + chip_cor + gain_cor
                    diff = mag_cor - mag0
                    _dmag_bad.append(diff)
                    _mag0_bad.append(mag0)
                    f.write("s %f %f %f %f %f 0\n" % (mag_cor, mag0, -9999,
                                                      s.u, s.v))
        f.close()

        d_mag_m = numpy.array(_dmag_m)
        d_mag_cat_m = numpy.array(_dmag_cat_m)
        d_mag_s = numpy.array(_dmag_s)
        d_mag_a = numpy.array(_dmag_a)
        d_mag_bad = numpy.array(_dmag_bad)
        d_mag_cat_bad = numpy.array(_dmag_cat_bad)
        mag0_m = numpy.array(_mag0_m)
        mag_cat_m = numpy.array(_mag_cat_m)
        mag0_s = numpy.array(_mag0_s)
        mag0_bad = numpy.array(_mag0_bad)
        mag_cat_bad = numpy.array(_mag_cat_bad)

        mag_std_m, mag_mean_m, mag_n_m  = self.clippedStd(d_mag_m, 3)
        mag_std_s, mag_mean_s, mag_n_s  = self.clippedStd(d_mag_s, 3)
        mag_std_a, mag_mean_a, mag_n_a  = self.clippedStd(d_mag_a, 3)
        mag_cat_std_m, mag_cat_mean_m, mag_cat_n_m  = self.clippedStd(d_mag_cat_m, 3)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
        plt.plot(mag0_bad, d_mag_bad, 'k,', markeredgewidth=0)
        plt.plot(mag_cat_m, d_mag_cat_m, 'c,', markeredgewidth=0)
        if self.sourceVec != None:
            plt.plot(mag0_s, d_mag_s, 'r,', markeredgewidth=0)
        plt.plot(mag0_m, d_mag_m, 'g,', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.xlim(15, 25)
        plt.ylim(-0.25, 0.25)
        plt.ylabel(r'$\Delta mag$ (mag)')

        bins = numpy.arange(-0.25, 0.25, 0.005) + 0.0025
        bins2 = numpy.arange(-0.25, 0.25, 0.05) + 0.025

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(d_mag_a, bins=bins, normed=False, orientation='horizontal', histtype='step')
        plt.hist(d_mag_m, bins=bins, normed=False, orientation='horizontal', histtype='step')
        if self.sourceVec != None:
            plt.hist(d_mag_s, bins=bins, normed=False, orientation='horizontal', histtype='step')
        plt.hist(d_mag_cat_m, bins=bins2, normed=False, orientation='horizontal', histtype='step')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (mag_std_a), rotation=270, transform=ax.transAxes, color='blue')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (mag_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.7, 0.90, r"$\sigma=$%5.3f" % (mag_cat_std_m), rotation=270, transform=ax.transAxes, color='cyan')
        y = mlab.normpdf(bins, mag_mean_m, mag_std_m)
        plt.plot(y*mag_n_m*0.005, bins, 'g:')
        if self.sourceVec != None:
            plt.text(0.3, 0.25, r"$\sigma=$%5.3f" % (mag_std_s), rotation=270, transform=ax.transAxes, color='red')
            y = mlab.normpdf(bins, mag_mean_s, mag_std_s)
            plt.plot(y*mag_n_s*0.005, bins, 'r:')
        y = mlab.normpdf(bins, mag_cat_mean_m, mag_cat_std_m)
        plt.plot(y*mag_cat_n_m*0.05, bins, 'c:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.25, 0.25)

        plt.savefig(os.path.join(self.outputDir, "MdM.png"), format='png')

    def plotPosDPos(self):
        _xi = []
        _eta = []
        _x = []
        _y = []
        for m in self.matchVec:
            if (m.good == True):
                _x.append((m.xi_fit - m.xi) * 3600)
                _y.append((m.eta_fit - m.eta) * 3600)
                _xi.append(m.xi * 3600)
                _eta.append(m.eta * 3600)
        if (self.sourceVec != None):
            for s in self.sourceVec:
                if (s.good == True):
                    _x.append((s.xi_fit - s.xi) * 3600)
                    _y.append((s.eta_fit - s.eta) * 3600)
                    _xi.append(s.xi * 3600)
                    _eta.append(s.eta * 3600)

        xi = numpy.array(_xi)
        eta = numpy.array(_eta)
        d_xi = numpy.array(_x)
        d_eta = numpy.array(_y)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot(2, 2, 1)
        plt.plot(xi, d_xi, ',', markeredgewidth=0)
        plt.xlabel(r'$\xi$ (arcsec)')
        plt.ylabel(r'$\Delta\xi$ (arcsec)')

        plt.subplot(2, 2, 3)
        plt.plot(xi, d_eta, ',', markeredgewidth=0)
        plt.xlabel(r'$\xi$ (arcsec)')
        plt.ylabel(r'$\Delta\eta$ (arcsec)')

        plt.subplot(2, 2, 2)
        plt.plot(eta, d_xi, ',', markeredgewidth=0)
        plt.xlabel(r'$\eta$ (arcsec)')
        plt.ylabel(r'$\Delta\xi$ (arcsec)')

        plt.subplot(2, 2, 4)
        plt.plot(eta, d_xi, ',', markeredgewidth=0)
        plt.xlabel(r'$\eta$ (arcsec)')
        plt.ylabel(r'$\Delta\eta$ (arcsec)')

        plt.savefig(os.path.join(self.outputDir, "PosDPos.png"), format='png')

    def plotResFlux(self):
        _dmag = []
        _iexp = []
        _ichip = []
        _r = []
        for m in self.matchVec:
            if (m.good == True and m.mag != -9999 and m.jstar != -1):
                mag = m.mag
                mag0 = m.mag0
                exp_cor = -2.5 * math.log10(self.fexp[m.iexp])
                chip_cor = -2.5 * math.log10(self.fchip[m.ichip])
                gain_cor = self.ffp.eval(m.u, m.v)
                mag_cor = mag + exp_cor + chip_cor + gain_cor
                diff = mag_cor - mag0
                _dmag.append(diff)
                _iexp.append(m.iexp)
                _ichip.append(m.ichip)

        d_mag = numpy.array(_dmag)
        iexp = numpy.array(_iexp)
        ichip = numpy.array(_ichip)

        mag_std  = self.clippedStd(d_mag, 3)[0]

        _r = []
        _dm = []
        for ccd in self.ccdSet.values():
            w = ccd.getAllPixels(True).getWidth()
            h = ccd.getAllPixels(True).getHeight()
            _x0 = ccd.getCenter().getPixels(ccd.getPixelSize())[0] + 0.5 * w
            _y0 = ccd.getCenter().getPixels(ccd.getPixelSize())[1] + 0.5 * h
            _r.append(math.sqrt(_x0*_x0 + _y0*_y0))
            _dm.append(-2.5 * math.log10(self.fchip[ccd.getId().getSerial()]))

        r = numpy.array(_r)
        dm = numpy.array(_dm)

        plt.clf()
        plt.rc('text', usetex=True)

        ax = plt.subplot(2, 2, 1)
        plt.hist(d_mag, bins=100, normed=True, histtype='step')
        plt.text(0.1, 0.7, r"$\sigma=$%7.5f" % (mag_std), transform=ax.transAxes)
        plt.xlabel(r'$\Delta mag$ (mag)')

        ax = plt.subplot(2, 2, 2)
        plt.plot(r, dm, 'o')
        plt.xlabel(r'Distance from center (pixel)')
        plt.ylabel(r'Offset in magnitude')

        ax = plt.subplot(2, 2, 3)
        plt.plot(iexp, d_mag, ',', markeredgewidth=0)
        plt.xlabel(r'Exposure ID')
        plt.ylabel(r'$\Delta mag$ (mag)')
        plt.xlim(iexp.min()-1, iexp.max()+1)
        plt.ylim(-0.2, 0.2)

        ax = plt.subplot(2, 2, 4)
        plt.plot(ichip, d_mag, ',', markeredgewidth=0)
        plt.xlabel(r'Chip ID')
        plt.ylabel(r'$\Delta mag$ (mag)')
        plt.xlim(ichip.min()-1, ichip.max()+1)
        plt.ylim(-0.2, 0.2)

        plt.savefig(os.path.join(self.outputDir, "ResFlux.png"), format='png')

    def plotDFlux2D(self):
        _dmag = []
        _u = []
        _v = []
        for m in self.matchVec:
            if (m.good == True and m.mag != -9999 and m.jstar != -1):
                mag = m.mag
                mag0 = m.mag0
                exp_cor = -2.5 * math.log10(self.fexp[m.iexp])
                chip_cor = -2.5 * math.log10(self.fchip[m.ichip])
                gain_cor = self.ffp.eval(m.u, m.v)
                mag_cor = mag + exp_cor + chip_cor + gain_cor
                diff = mag_cor - mag0
                _dmag.append(diff)
                _u.append(m.u)
                _v.append(m.v)

        d_mag = numpy.array(_dmag)
        u = numpy.array(_u)
        v = numpy.array(_v)

        s = numpy.absolute(d_mag) * 10

        u1 = [u[i] for i in range(len(d_mag)) if d_mag[i] > 0]
        v1 = [v[i] for i in range(len(d_mag)) if d_mag[i] > 0]
        s1 = [math.fabs(d_mag[i])*20 for i in range(len(d_mag)) if d_mag[i] > 0]
        u2 = [u[i] for i in range(len(d_mag)) if d_mag[i] < 0]
        v2 = [v[i] for i in range(len(d_mag)) if d_mag[i] < 0]
        s2 = [math.fabs(d_mag[i])*20 for i in range(len(d_mag)) if d_mag[i] < 0]

        plt.clf()
        plt.rc('text', usetex=True)

        plt.scatter(u1, v1, s1, color='blue')
        plt.scatter(u2, v2, s2, color='red')
        plt.axes().set_aspect('equal')

        plt.savefig(os.path.join(self.outputDir, "DFlux2D.png"), format='png')

    def outputDiag(self):
        self.log.info("Output Diagnostic Figures...")

        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        f = open(os.path.join(self.outputDir, "coeffs.dat"), "wt")
        for iexp in self.coeffSet.keys():
            c = self.coeffSet[iexp]
            f.write("%ld %12.5e %12.5e\n" % (iexp, c.A,  c.D));
            f.write("%ld %12.5f %12.5f\n" % (iexp, c.x0, c.y0));
            for k in range(c.getNcoeff()):
                f.write("%ld %15.8e %15.8e %15.8e %15.8e\n" % (iexp, c.get_a(k), c.get_b(k), c.get_ap(k), c.get_bp(k)));
            f.write("%5.3f\n" % (-2.5*math.log10(self.fexp[iexp])))
        f.close()

        f = open(os.path.join(self.outputDir, "ccd.dat"), "wt")
        for ichip in self.ccdSet.keys():
            ccd = self.ccdSet[ichip]
            center = ccd.getCenter().getPixels(ccd.getPixelSize())
            orient = ccd.getOrientation()
            f.write("%3ld %10.3f %10.3f %10.7f %5.3f\n" % (ichip, center[0], center[1], orient.getYaw(), self.fchip[ichip]));
        f.close()

        for iexp in self.coeffSet.keys():
            self.plotJCont(iexp)
            self.plotFCorCont(iexp)
            self.plotResPosArrow2D(iexp)

        self.plotResPosScatter()
        self.plotMdM()
        self.plotPosDPos()
        self.plotResFlux()
        self.plotDFlux2D()

    def mosaic(self, butler, frameIds, ccdIds, ct=None, debug=False, verbose=False):

        self.log.info(str(self.config))

        ccdSet = self.readCcd(butler.mapper.camera, ccdIds)

        if debug:
            for ccd in ccdSet.values():
                self.log.info(str(ccd.getId().getSerial())+" "+
                              str(ccd.getCenter().getPixels(ccd.getPixelSize()))+" "+
                              str(ccd.getOrientation().getYaw()))

        wcsDic = self.readWcs(butler, frameIds, ccdSet)

        self.removeNonExistCcd(butler, ccdSet, wcsDic)

        if debug:
            for iexp, wcs in wcsDic.iteritems():
                self.log.info(str(iexp)+" "+str(wcs.getPixelOrigin())+" "+
                              str(wcs.getSkyOrigin().getPosition(afwGeom.degrees)))

        sourceSet, matchList = self.readCatalog(butler, wcsDic.keys(), ccdSet.keys(), ct)

        wcsDic, sourceSet, matchList = self.checkInputs(wcsDic, sourceSet, matchList)
        self.log.info("frameIds : "+str(wcsDic.keys()))
        self.log.info("ccdIds : "+str(ccdSet.keys()))

        d_lim = afwGeom.Angle(self.config.radXMatch, afwGeom.arcseconds)
        nbrightest = self.config.nBrightest
        if debug:
            self.log.info("d_lim : %f" % d_lim)
            self.log.info("nbrightest : %d" % nbrightest)

        allMat, allSource =self.mergeCatalog(sourceSet, matchList, ccdSet, d_lim, nbrightest)

        nmatch  = allMat.size()
        nsource = allSource.size()
        matchVec  = measMosaic.obsVecFromSourceGroup(allMat,    wcsDic, ccdSet)
        sourceVec = measMosaic.obsVecFromSourceGroup(allSource, wcsDic, ccdSet)

        self.log.info("Solve mosaic ...")
        order = self.config.fittingOrder
        internal = self.config.internalFitting
        solveCcd = self.config.solveCcd
        allowRotation = self.config.allowRotation
        fluxFitOrder = self.config.fluxFitOrder
        chebyshev = self.config.chebyshev
        absolute = self.config.fluxFitAbsolute
        catRMS = self.config.catRMS

        if not internal:
            sourceVec = None

        if debug:
            self.log.info("order : %d" % ffp.order)
            self.log.info("internal : %r" % internal)
            self.log.info("solveCcd : %r " % solveCcd)
            self.log.info("allowRotation : %r" % allowRotation)

        ffp = measMosaic.FluxFitParams(fluxFitOrder, absolute, chebyshev)
        u_max, v_max = self.getExtent(matchVec)
        ffp.u_max = (math.floor(u_max / 10.) + 1) * 10
        ffp.v_max = (math.floor(v_max / 10.) + 1) * 10

        fexp = measMosaic.map_int_float()
        fchip = measMosaic.map_int_float()

        if internal:
            coeffSet = measMosaic.solveMosaic_CCD(order, nmatch, nsource,
                                                  matchVec, sourceVec,
                                                  wcsDic, ccdSet, ffp, fexp, fchip,
                                                  solveCcd, allowRotation, verbose, catRMS)
        else:
            coeffSet = measMosaic.solveMosaic_CCD_shot(order, nmatch, matchVec, 
                                                       wcsDic, ccdSet, ffp, fexp, fchip,
                                                       solveCcd, allowRotation, verbose, catRMS)

        self.butler = butler
        self.outputDir = self.config.outputDir
        self.matchVec = matchVec
        self.sourceVec = sourceVec
        self.wcsDic = wcsDic
        self.ccdSet = ccdSet
        self.coeffSet = coeffSet
        self.ffp = ffp
        self.fexp = fexp
        self.fchip = fchip

        self.writeNewWcs()
        self.writeFcr()

        if self.config.outputDiag:
            self.outputDiag()

        return wcsDic.keys()

    def run(self, camera, butler, dataRefList, debug, verbose=False):

        frameIds = list()
        ccdIds = list()
        filters = list()
        fields = list()
        for dataRef in dataRefList:
            if not dataRef.dataId['visit'] in frameIds:
                frameIds.append(dataRef.dataId['visit'])
            if not dataRef.dataId['ccd'] in ccdIds:
                ccdIds.append(dataRef.dataId['ccd'])
            if not dataRef.dataId['filter'] in filters:
                filters.append(dataRef.dataId['filter'])
            if not dataRef.dataId['field'] in fields:
                fields.append(dataRef.dataId['field'])

        if len(filters) != 1:
            self.log.fatal("There are %d filters in input frames" % len(filters))
            return None

        if camera == 'suprimecam':
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("Hamamatsu")
            ct = Colorterm.getColorterm(butler.mapper.filters[filters[0]])
        elif camera == 'suprimecam-mit':
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("MIT")
            ct = Colorterm.getColorterm(butler.mapper.filters[filters[0]])
        else:
            ct = None

        return self.mosaic(butler, frameIds, ccdIds, ct, debug, verbose)

    def getAllForCcdNew(self, butler, astrom, frame, ccd, ct=None):

        data = {'visit': frame, 'ccd': ccd}

        try:
            if not butler.datasetExists('src', data):
                raise RuntimeError("no data for src %s" % (data))
            if not butler.datasetExists('calexp_md', data):
                raise RuntimeError("no data for calexp_md %s" % (data))

            wcs_md = butler.get('wcs_md', data)
            wcs = afwImage.makeWcs(wcs_md)

            fcr_md = butler.get('fcr_md', data)
            ffp = measMosaic.FluxFitParams(fcr_md)
            fluxmag0 = fcr_md.get('FLUXMAG0')

            sources = butler.get('src', data)
            if False:
                matches = measAstrom.readMatches(butler, data)
            else:
                icSrces = butler.get('icSrc', data)
                packedMatches = butler.get('icMatch', data)
                matches = astrom.joinMatchListWithCatalog(packedMatches, icSrces, True)
                if ct != None:
                    if matches[0].first != None:
                        refSchema = matches[0].first.schema
                    else:
                        refSchema = matches[1].first.schema
                    key_p = refSchema.find(ct.primary).key
                    key_s = refSchema.find(ct.secondary).key
                    key_f = refSchema.find("flux").key
                    for m in matches:
                        if m.first != None:
                            refFlux1 = m.first.get(key_p)
                            refFlux2 = m.first.get(key_s)
                            refMag1 = -2.5*math.log10(refFlux1)
                            refMag2 = -2.5*math.log10(refFlux2)
                            refMag = ct.transformMags(ct.primary, refMag1, refMag2)
                            refFlux = math.pow(10.0, -0.4*refMag)
                            if refFlux == refFlux:
                                m.first.set(key_f, refFlux)
                            else:
                                m.first = None

            sources = self.selectStars(sources)
            selMatches = self.selectStars(matches)
            if len(selMatches) < 10:
                matches = self.selectStars(matches, True)
            else:
                matches = selMatches
        except Exception, e:
            print "Failed to read: %s" % (e)
            return None, None, None, None, None
    
        return sources, matches, wcs, ffp, fluxmag0

    def check(self, butler, frameIds, ccdIds, ct=None, debug=False, verbose=False):

        self.log.info(str(self.config))

        astrom = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
        for frameId in frameIds:
            for ccdId in ccdIds:
                sources, matches, wcs, ffp, fluxmag0 = self.getAllForCcdNew(butler, astrom, frameId, ccdId, ct)
                if matches != None and len(matches) != 0:
                    refSchema = matches[0].first.schema
                    key_p = refSchema.find(ct.primary).key
                    key_s = refSchema.find(ct.secondary).key
                    key_f = refSchema.find("flux").key
                    key_e = refSchema.find("flux.err").key
                    for m in matches:
                        cat = m.first
                        src = m.second
                        if cat != None and src != None:
                            src_sky = wcs.pixelToSky(src.getX(), src.getY())
                            src.setRa(src_sky[0])
                            src.setDec(src_sky[1])
                            mag_cat = -2.5*math.log10(cat.get(key_f))
                            mag_err = 2.5/math.log(10)*cat.get(key_e)/cat.get(key_f)
                            mag_p = -2.5*math.log10(cat.get(key_p))
                            mag_s = -2.5*math.log10(cat.get(key_s))
                            #mag_src = -2.5*math.log10(src.getPsfFlux()/fluxmag0)
                            try:
                                mag_src = -2.5*math.log10(src.getApFlux()/fluxmag0)
                                print frameId, ccdId, mag_cat, \
                                    (src.getRa()-cat.getRa()).asArcseconds(), \
                                    (src.getDec()-cat.getDec()).asArcseconds(), \
                                    mag_src-mag_cat+ffp.eval(src.getX(), src.getY()), \
                                    mag_p, mag_s, mag_err, ffp.eval(src.getX(), src.getY())
                            except Exception, e:
                                pass

        return None

    def makeDiffPos(self, allMat, allSource, wcsDic):
        dx_m = list()
        dy_m = list()
        dx_s = list()
        dy_s = list()
        for ss in allMat:
            ra_cat = ss[0].getRa().asDegrees()
            dec_cat = ss[0].getDec().asDegrees()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                x = ss[j].getX()
                y = ss[j].getY()
                sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                ra = sky[0]
                dec = sky[1]
                dx_m.append((ra - ra_cat) * 3600)
                dy_m.append((dec - dec_cat) * 3600)
            if len(ss) > 2:
                n = 0
                ra_cat = 0.0
                dec_cat = 0.0
                ra_source = list()
                dec_source = list()
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                    ra = sky[0]
                    dec = sky[1]
                    ra_source.append(ra)
                    dec_source.append(dec)
                    n += 1
                    ra_cat += ra
                    dec_cat += dec
                ra_cat /= n
                dec_cat /= n
                for ra, dec in zip(ra_source, dec_source):
                    dx_s.append((ra - ra_cat) * 3600)
                    dy_s.append((dec - dec_cat) * 3600)
        dx_m = numpy.array(dx_m)
        dy_m = numpy.array(dy_m)

        for ss in allSource:
            n = 0
            ra_cat = 0.0
            dec_cat = 0.0
            ra_source = list()
            dec_source = list()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                x = ss[j].getX()
                y = ss[j].getY()
                sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                ra = sky[0]
                dec = sky[1]
                ra_source.append(ra)
                dec_source.append(dec)
                n += 1
                ra_cat += ra
                dec_cat += dec
            ra_cat /= n
            dec_cat /= n
            for ra, dec in zip(ra_source, dec_source):
                dx_s.append((ra - ra_cat) * 3600)
                dy_s.append((dec - dec_cat) * 3600)
        dx_s = numpy.array(dx_s)
        dy_s = numpy.array(dy_s)

        return dx_m, dy_m, dx_s, dy_s

    def makeDiffFlux(self, allMat, allSource, calibDic, ffpDic, mag_lim = 9999.0):

        mag0_m = list()
        mag_m  = list()
        mcor_m = list()
        mag0_s = list()
        mag_s  = list()
        mcor_s = list()
        for ss in allMat:
            mag_cat = -2.5*math.log10(ss[0].getFlux())
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                if (ss[j].getFlux() > 0):
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag0_m.append(mag_cat)
                    mag_m.append(mag)
                    mcor_m.append(mcor)
            if len(ss) > 2:
                Sx = 0.0
                S  = 0.0
                mag_source = list()
                mcor_source = list()
                iexp_source = list()
                ichip_source = list()
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag_source.append(mag)
                    mcor_source.append(mcor)
                    iexp_source.append(iexp)
                    ichip_source.append(ichip)
                    Sx += (mag+mcor) / (err*err)
                    S  += 1. / (err*err)
                mag_cat = Sx / S
                if mag_cat < mag_lim:
                    for mag, mcor, iexp, ichip in zip(mag_source, mcor_source, iexp_source, ichip_source):
                        mag0_s.append(mag_cat)
                        mag_s.append(mag)
                        mcor_s.append(mcor)
        mag0_m = numpy.array(mag0_m)
        mag_m  = numpy.array(mag_m)
        mcor_m = numpy.array(mcor_m)
        dm_m = mag_m + mcor_m - mag0_m

        for ss in allSource:
            Sx = 0.0
            S  = 0.0
            mag_source = list()
            mcor_source = list()
            iexp_source = list()
            ichip_source = list()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                #print iexp, ichip, calibDic[iexp][ichip].getFluxMag0()[0], ss[j].getFlux()
                if calibDic[iexp][ichip].getFluxMag0()[0] > 0:
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag_source.append(mag)
                    mcor_source.append(mcor)
                    iexp_source.append(iexp)
                    ichip_source.append(ichip)
                    Sx += (mag+mcor) / (err*err)
                    S  += 1. / (err*err)
            if S > 0:
                mag_cat = Sx / S
            else:
                mag_cat = 99999.
            if mag_cat < mag_lim:
                for mag, mcor, iexp, ichip in zip(mag_source, mcor_source, iexp_source, ichip_source):
                    mag0_s.append(mag_cat)
                    mag_s.append(mag)
                    mcor_s.append(mcor)
        mag0_s = numpy.array(mag0_s)
        mag_s  = numpy.array(mag_s)
        mcor_s = numpy.array(mcor_s)
        dm_s = mag_s + mcor_s - mag0_s

        return mag0_m, dm_m, mag0_s, dm_s

    def makeFluxStat(self, allMat, allSource, calibDic, ffpDic, mag_lim = 9999.0):

        x = list()
        y = list()
        ra = list()
        dec = list()
        id = list()
        for ss in allMat:
            Sxx = 0.0
            Sx  = 0.0
            S   = 0.0
            Sr  = 0.0
            Sd  = 0.0
            if len(ss) > 2:
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    mcor = ffpDic[iexp][ichip].eval(ss[j].getX(), ss[j].getY())
                    Sxx += (mag+mcor)*(mag+mcor) / (err*err)
                    Sx  += (mag+mcor) / (err*err)
                    S   += 1. / (err*err)
                    Sr += ss[j].getRa().asDegrees()
                    Sd += ss[j].getDec().asDegrees()
                avg = Sx / S
                sig = math.sqrt(Sxx/S - avg*avg)
                x.append(avg)
                y.append(sig)
                ra.append(Sr / (len(ss)-1))
                dec.append(Sd / (len(ss)-1))

        for ss in allSource:
            Sxx = 0.0
            Sx = 0.0
            S  = 0.0
            Sr  = 0.0
            Sd  = 0.0
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                #print iexp, ichip, calibDic[iexp][ichip].getFluxMag0()[0], ss[j].getFlux()
                if calibDic[iexp][ichip].getFluxMag0()[0] > 0:
                    mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    mcor = ffpDic[iexp][ichip].eval(ss[j].getX(), ss[j].getY())
                    Sxx += (mag+mcor)*(mag+mcor) / (err*err)
                    Sx  += (mag+mcor) / (err*err)
                    S   += 1. / (err*err)
                    Sr += ss[j].getRa().asDegrees()
                    Sd += ss[j].getDec().asDegrees()
            if S > 0:
                try:
                    avg = Sx / S
                    sig = math.sqrt(Sxx/S - avg*avg)
                    x.append(avg)
                    y.append(sig)
                    ra.append(Sr / (len(ss)-1))
                    dec.append(Sd / (len(ss)-1))
                except Exception, e:
                    #print Sxx, S, avg, Sxx/S - avg*avg, len(ss)-1
                    pass

        if True:
            plt.clf()
            plt.plot(x, y, ',', markeredgewidth=0)
            plt.xlim(15, 25)
            plt.ylim(0.0, 0.20)
            plt.plot([15, 25], [0.01, 0.01], 'k--')
            plt.xlabel('mag (avg)')
            plt.ylabel('RMS')
            #plt.title('r-band')
            plt.savefig('fluxMean.png')
        else:
            for r, d, m, dm in zip(ra, dec, x, y):
                print '%9.5f %9.5f %7.4f %7.4f' % (r, d, m ,dm)

    def plotPos(self, dx_m, dy_m, dx_s, dy_s):

        x_std_m, x_mean_m, x_n_m = self.clippedStd(dx_m, 3)
        y_std_m, y_mean_m, y_n_m = self.clippedStd(dy_m, 3)
        x_std_s, x_mean_s, x_n_s = self.clippedStd(dx_s, 3)
        y_std_s, y_mean_s, y_n_s = self.clippedStd(dy_s, 3)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6), (1,0), colspan=4, rowspan=4)
        plt.plot(dx_m, dy_m, 'g,', markeredgecolor='green')
        plt.plot(dx_s, dy_s, 'r,', markeredgecolor='red')
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        plt.xlabel(r'$\Delta\alpha$ (arcsec)')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')

        bins = numpy.arange(-0.5, 0.5, 0.01) + 0.005

        ax = plt.subplot2grid((5,6),(0,0), colspan=4)
        plt.hist([dx_m, dx_s], bins=bins, normed=False, histtype='step', color=['green', 'red'])
        plt.text(0.75, 0.7, r"$\sigma=$%5.3f" % (x_std_m), transform=ax.transAxes, color='green')
        plt.text(0.75, 0.5, r"$\sigma=$%5.3f" % (x_std_s), transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, x_mean_m, x_std_m)
        plt.plot(bins, gauss*x_n_m*0.01, 'g:')
        gauss = mlab.normpdf(bins, x_mean_s, x_std_s)
        plt.plot(bins, gauss*x_n_s*0.01, 'r:')
        plt.xlim(-0.5, 0.5)

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dy_m, bins=bins, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.hist(dy_s, bins=bins, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (y_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (y_std_s), rotation=270, transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, y_mean_m, y_std_m)
        plt.plot(gauss*y_n_m*0.01, bins, 'g:')
        gauss = mlab.normpdf(bins, y_mean_s, y_std_s)
        plt.plot(gauss*y_n_s*0.01, bins, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.5, 0.5)

        plt.savefig("posScatter.png")

    def plotFlux(self, m0_m, dm_m, m0_s, dm_s):

        mag_std_m, mag_mean_m, mag_n_m = self.clippedStd(dm_m, 3)
        mag_std_s, mag_mean_s, mag_n_s = self.clippedStd(dm_s, 3)

        bins = numpy.arange(-0.25, 0.25, 0.005) + 0.0025

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
        plt.plot(m0_m, dm_m, 'g,', markeredgecolor='green')
        plt.plot(m0_s, dm_s, 'r,', markeredgecolor='red')
        plt.plot([15,25], [0,0], 'k--')
        #binnedScat(m0_s, dm_s)
        plt.xlim(15, 25)
        plt.ylim(-0.25, 0.25)
        plt.plot([15, 25], [-0.01, -0.01], 'k--')
        plt.plot([15, 25], [+0.01, +0.01], 'k--')
        plt.xlabel(r'$m_{cat}$ (mag)')
        plt.ylabel(r'$\Delta m$ (mag)')

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dm_s, bins=bins, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.hist(dm_m, bins=bins, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (mag_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (mag_std_s), rotation=270, transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, mag_mean_m, mag_std_m)
        plt.plot(gauss*mag_n_m*0.005, bins, 'g:')
        gauss = mlab.normpdf(bins, mag_mean_s, mag_std_s)
        plt.plot(gauss*mag_n_s*0.005, bins, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.25, 0.25)

        plt.savefig("fluxScatter.png")

    def plotPosAsMag(self, m0_s, dx_s, dy_s):

        plt.subplot2grid((2,5),(0,0), colspan=4)
        plt.plot(m0_s, dx_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\alpha$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        mlim = 24
        ax = plt.subplot2grid((2,5),(0,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dx_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dx_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.subplot2grid((2,5),(1,0), colspan=4)
        plt.plot(m0_s, dy_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        ax = plt.subplot2grid((2,5),(1,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dy_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dy_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.savefig('posAsMag.png')

    def check2(self, butler, frameIds, ccdIds, ct=None, debug=False, verbose=False):
        self.log.info(str(self.config))

        ccdSet = self.readCcd(butler.mapper.camera, ccdIds)
        sourceSet = measMosaic.SourceGroup()
        matchList = measMosaic.SourceMatchGroup()
        wcsDic = dict()
        calibDic = dict()
        ffpDic = dict()
        astrom = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
        for frameId in frameIds:
            data = {'visit': frameId, 'ccd': 0}
            if not butler.datasetExists('src', data):
                self.log.info(str(data)+" is not exist")
                continue
            ss = []
            ml = []
            if not wcsDic.has_key(frameId):
                wcsDic[frameId] = dict()
            if not calibDic.has_key(frameId):
                calibDic[frameId] = dict()
            if not ffpDic.has_key(frameId):
                ffpDic[frameId] = dict()
            for ccdId in ccdIds:

                #sources, matches, wcs, calib, ffp = getAllForCcd(butler, astrom, frameId, ccdId)

                data = {'visit': frameId, 'ccd': ccdId}
                self.log.info(str(data))
                try:
                    if not butler.datasetExists('src', data):
                        raise RuntimeError("no data for src %s" % (data))
                    #if not butler.datasetExists('calexp_md', data):
                    #    raise RuntimeError("no data for calexp_md %s" % (data))
                    #md = butler.get('calexp_md', data)
                    md = butler.get('wcs_md', data)
                    wcs = afwImage.makeWcs(md)
                    #wcs = afwImage.makeWcs(wcs.getFitsMetadata())
                    calib = afwImage.Calib(md)

                    md = butler.get('fcr_md', data)
                    ffp = measMosaic.FluxFitParams(md)

                    sources = butler.get('src', data)
                    if False:
                        matches = measAstrom.readMatches(butler, data)
                    else:
                        icSrces = butler.get('icSrc', data)
                        packedMatches = butler.get('icMatch', data)
                        matches = astrom.joinMatchListWithCatalog(packedMatches, icSrces, True)

                    if ct != None:
                        if matches[0].first != None:
                            refSchema = matches[0].first.schema
                        else:
                            refSchema = matches[1].first.schema
                        key_p = refSchema.find(ct.primary).key
                        key_s = refSchema.find(ct.secondary).key
                        key_f = refSchema.find("flux").key
                        for m in matches:
                            if m.first != None:
                                refFlux1 = m.first.get(key_p)
                                refFlux2 = m.first.get(key_s)
                                refMag1 = -2.5*math.log10(refFlux1)
                                refMag2 = -2.5*math.log10(refFlux2)
                                refMag = ct.transformMags(ct.primary, refMag1, refMag2)
                                refFlux = math.pow(10.0, -0.4*refMag)
                                if refFlux == refFlux:
                                    m.first.set(key_f, refFlux)
                                else:
                                    m.first = None

                    sources = self.selectStars(sources)
                    matches = self.selectStars(matches, True)
                except Exception, e:
                    print "Failed to read: %s" % (e)

                if sources != None:
                    for s in sources:
                        if numpy.isfinite(s.getRa().asDegrees()): # get rid of NaN
                            src = measMosaic.Source(s)
                            src.setExp(frameId)
                            src.setChip(ccdId)
                            ss.append(src)
                    for m in matches:
                        match = measMosaic.SourceMatch(measMosaic.Source(m.first, wcs),
                                                         measMosaic.Source(m.second))
                        match.second.setExp(frameId)
                        match.second.setChip(ccdId)
                        ml.append(match)
                wcsDic[frameId][ccdId] = wcs
                calibDic[frameId][ccdId] = calib
                ffpDic[frameId][ccdId] = ffp
            sourceSet.push_back(ss)
            matchList.push_back(ml)

        d_lim = afwGeom.Angle(self.config.radXMatch, afwGeom.arcseconds)
        nbrightest = self.config.nBrightest
        allMat, allSource = self.mergeCatalog(sourceSet, matchList, ccdSet, d_lim, nbrightest)

        m0_m, dm_m, m0_s, dm_s = self.makeDiffFlux(allMat, allSource, calibDic, ffpDic)
        print len(m0_m), len(dm_m), len(m0_s), len(dm_s)
        self.plotFlux(m0_m, dm_m, m0_s, dm_s)
        self.makeFluxStat(allMat, allSource, calibDic, ffpDic)
        dx_m, dy_m, dx_s, dy_s = self.makeDiffPos(allMat, allSource, wcsDic)
        print len(dx_m), len(dy_m), len(dx_s), len(dy_s)
        self.plotPos(dx_m, dy_m, dx_s, dy_s)
        self.plotPosAsMag(m0_s, dx_s, dy_s)

    def run_check(self, camera, butler, dataRefList, debug, verbose=False):

        frameIds = list()
        ccdIds = list()
        filters = list()
        for dataRef in dataRefList:
             if not dataRef.dataId['visit'] in frameIds:
                 frameIds.append(dataRef.dataId['visit'])
             if not dataRef.dataId['ccd'] in ccdIds:
                 ccdIds.append(dataRef.dataId['ccd'])
             if not dataRef.dataId['filter'] in filters:
                 filters.append(dataRef.dataId['filter'])

        if len(filters) != 1:
            self.log.fatal("There are %d filters in input frames" % len(filters))
            return None

        if camera == 'suprimecam':
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("Hamamatsu")
            ct = Colorterm.getColorterm(butler.mapper.filters[filters[0]])
        elif camera == 'suprimecam-mit':
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("MIT")
            ct = Colorterm.getColorterm(butler.mapper.filters[filters[0]])
        else:
            ct = None

        #return self.check(butler, frameIds, ccdIds, ct, debug, verbose)
        return self.check2(butler, frameIds, ccdIds, ct, debug, verbose)
