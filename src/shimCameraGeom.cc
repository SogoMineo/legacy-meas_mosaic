#include "lsst/meas/mosaic/shimCameraGeom.h"
#include "lsst/pex/exceptions.h"

namespace lsst {
namespace meas {
namespace mosaic {

int getNQuarter(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getOrientation().getNQuarter();
}

lsst::geom::Angle getYaw(CONST_PTR(afw::cameraGeom::Detector) det) {
    lsst::geom::Angle deg = det->getOrientation().getYaw();
    int nQuarter = det->getOrientation().getNQuarter();
    if (nQuarter%4 != 0) {
        deg = det->getOrientation().getYaw() - nQuarter*90.0*lsst::geom::degrees;
    }
    if (fabs(deg.asDegrees()) >= 90.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          (boost::format("Mismatch between yaw (%f deg) and nQuarter (%d) for detector %d:"
                           " abs(yaw - 90*nQuarter) = %f is > 90 deg")
                          % det->getOrientation().getYaw().asDegrees()
                          % getNQuarter(det)
                          % det->getSerial()
                          % fabs(deg.asDegrees())).str());
    }
    return deg;
}

lsst::geom::LinearTransform makeScalingMmToPx(lsst::geom::Extent2D const pSize) {
    return lsst::geom::LinearTransform::makeScaling(1.0/pSize.getX(), 1.0/pSize.getY());
}

lsst::geom::Point2D getCenterInFpPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->getCenter(afw::cameraGeom::FOCAL_PLANE));
}

lsst::geom::Point2D getCenterInDetectorPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto center = det->getCenter(afw::cameraGeom::PIXELS);
    if ((getNQuarter(det)%2) != 0) {
        return lsst::geom::Point2D(center.getY(), center.getX());
    } else {
        return center;
    }
}

int getWidth(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getBBox().getWidth();
}

int getHeight(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getBBox().getHeight();
}

lsst::geom::Point2D detPxToFpPx(CONST_PTR(afw::cameraGeom::Detector) det, lsst::geom::Point2D const detPt) {
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->transform(detPt, afw::cameraGeom::PIXELS, afw::cameraGeom::FOCAL_PLANE));
}

lsst::geom::Point2D detPxToFpPxRot(CONST_PTR(afw::cameraGeom::Detector) det, lsst::geom::Point2D const detPt) {
    double cosYaw = std::cos(getYaw(det));
    double sinYaw = std::sin(getYaw(det));
    // Center in detector and focal plane pixels
    lsst::geom::Point2D centerDet = getCenterInDetectorPixels(det);
    lsst::geom::Point2D centerFp = getCenterInFpPixels(det);

    lsst::geom::Extent2D offset = lsst::geom::Extent2D(cosYaw * detPt.getX() - sinYaw * detPt.getY(),
                                                     sinYaw * detPt.getX() + cosYaw * detPt.getY());
    offset -= lsst::geom::Extent2D(centerDet);
    return centerFp + offset;
}

lsst::geom::Point2D computeX0Y0(CONST_PTR(afw::cameraGeom::Detector) det, double x0, double y0) {
    lsst::geom::Point2D newXY0;

    double cosYaw = std::cos(getYaw(det));
    double sinYaw = std::sin(getYaw(det));

    // Offset between center in focal plane and detector pixels
    lsst::geom::Extent2D off = getCenterInFpPixels(det) - getCenterInDetectorPixels(det);

    newXY0[0] =  (off[0] + x0)*cosYaw + (off[1] + y0)*sinYaw;
    newXY0[1] = -(off[0] + x0)*sinYaw + (off[1] + y0)*cosYaw;

    return newXY0;
}

}}} // lsst::meas::mosaic
