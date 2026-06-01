/*
 *  SPDX-FileCopyrightText: 2026 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "LcmsPredefinedPipelineFunctions.h"
#include <KoColorTransferFunctions.h>
#include <KoColorProfile.h>
#include <QMatrix4x4>
#include <QVector4D>

/**
 * @brief bradfordMatrix
 * From http://brucelindbloom.com/index.html?Eqn_ChromAdapt.html
 * @param src -- source white point XYZ
 * @param dst -- destionation white point XYZ
 * @return chromatic adaptation matrix.
 */
static QMatrix4x4 bradfordMatrix(QVector4D src, QVector4D dst) {

    QVector4D srcLMS(src[0], src[1], src[2], 0.0);
    QVector4D dstLMS(dst[0], dst[1], dst[2], 0.0);

    const QMatrix4x4 bradfordMatrix(
        0.8951000, 0.2664000, -0.1614000, 0.0,
        -0.7502000, 1.7135000, 0.0367000, 0.0,
        0.0389000, -0.0685000, 1.0296000, 0.0,
        0.0, 0.0, 0.0, 1.0
        );
    const QMatrix4x4 bradfordMatrixInverted (
        0.9869929, -0.1470543, 0.1599627, 0.0,
        0.4323053, 0.5183603, 0.0492912, 0.0,
        -0.0085287, 0.0400428, 0.9684867, 0.0,
        0.0, 0.0, 0.0, 1.0
    );

    srcLMS = bradfordMatrix * srcLMS;
    dstLMS = bradfordMatrix * dstLMS;

    const QMatrix4x4 adaptationMatrix (
        dstLMS.x()/srcLMS.x(), 0, 0, 0,
        0, dstLMS.y()/srcLMS.y(), 0, 0,
        0, 0, dstLMS.z()/srcLMS.z(), 0,
        0.0, 0.0, 0.0, 1.0
    );


    return bradfordMatrixInverted * adaptationMatrix * bradfordMatrix;
}
/**
 * @brief rgbMatrix
 * Generates a proper rgb matrix from a primaries type.
 * @param primaries -- default primaries to use.
 * @param factor -- scaling factor.
 * @return transfrom representing a XYZ to RGB conversion.
 */
static QMatrix4x4 rgbMatrix(ColorPrimaries primaries, double factor = 1.0, double *luma = nullptr) {
    const double mul = 1.0/factor;
    // Calculation: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    // Calculates XYZ with Y = 1.0
    QVector<double> col;
    KoColorProfile::colorantsForType(primaries, col);
    QMatrix4x4 initialXYZMatrix(
        col[2]/col[3], col[4]/col[5], col[6]/col[7], 0.0,
        1.0, 1.0, 1.0, 0.0,
        (1.0-col[2]-col[3])/col[3], (1.0-col[4]-col[5])/col[5], (1.0-col[6]-col[7])/col[7], 0.0,
        0.0, 0.0, 0.0, 1.0
        );

    // Calculate scaling factor with white point.
    QVector4D wp(col[0]/col[1], 1.0, (1.0-col[0]-col[1])/col[1], 0.0);

    QVector4D S = initialXYZMatrix.inverted() * wp;
    if (luma) {
        luma[0] = S.x();
        luma[1] = S.y();
        luma[2] = S.z();
    }

    // Apply scaling factor to normalize the Y and get the matrix.
    QMatrix4x4 finalXYZMatrix(
        initialXYZMatrix(0, 0) * S.x(), initialXYZMatrix(0, 1) * S.y(), initialXYZMatrix(0, 2) * S.z(), 0.0,
        initialXYZMatrix(1, 0) * S.x(), initialXYZMatrix(1, 1) * S.y(), initialXYZMatrix(1, 2) * S.z(), 0.0,
        initialXYZMatrix(2, 0) * S.x(), initialXYZMatrix(2, 1) * S.y(), initialXYZMatrix(2, 2) * S.z(), 0.0,
        0.0, 0.0, 0.0, 1.0
        );
    // Apply bradford white point adaptation.
    QVector4D wpD50(0.9642, 1.0, 0.8249, 0.0);
    finalXYZMatrix = bradfordMatrix(wp, wpD50) * finalXYZMatrix;

    finalXYZMatrix.scale(mul);
    return finalXYZMatrix;
}

struct perceptualDummyHelper {
    bool toXYZ{false};
    int channels {3};
    double *luma{nullptr};
    };

/**
 * @brief samplePQDummyClut
 * In this sampler, we convert the PQ data to linear and then to HLG. The benefit of HLG
 * is that it sort of allows for tonemapping.
 */
cmsInt32Number samplePQDummyClut(const cmsUInt16Number In[], cmsUInt16Number Out[], void *Cargo) {
    struct perceptualDummyHelper *helper = (struct perceptualDummyHelper *) Cargo;
    const float pqScale = 125.0; /// Important: this normalizes the pq signal.
    const float nominalPeak = 10000.0/392.0;/// See bt. 2390 pg. 52.
    double coeff[] = {0.2126, 0.7152, 0.0722};
    if (helper->luma) {
        coeff[0] = helper->luma[0];
        coeff[1] = helper->luma[1];
        coeff[2] = helper->luma[2];
    }
    float lin[3];
    for (int i = 0; i < helper->channels; i++) {
        const float val = float(In[i])/65535.0;
        if (helper->toXYZ) {
            lin[i] = (removeSmpte2048Curve(val) / pqScale) * nominalPeak;
        } else {
            lin[i] = removeHLGCurve(val);
        }
    }
    for (int i = 0; i < helper->channels; i++) {
        cmsUInt16Number finalResult = In[i];
        if (helper->toXYZ) {
            finalResult = qBound(0, qRound( applyHLGCurve(lin[i]) * 65535 ), 65535);
        } else {
            finalResult = qBound(0, qRound( applySmpte2048Curve((lin[i] / nominalPeak) * pqScale) * 65535 ), 65535);
        }
        Out[i] = finalResult;
    }
    return true;
};

bool LcmsPredefinedPipelineFunctions::setPerceptualQuantizerAToBDummyPipeline(cmsHPROFILE iccProfile, cmsTagSignature tag, ColorPrimaries primaries)
{
    // From profile space to XYZ.
    const double factor = 1.0;
    const int lutSize = 4;

    double luma[3];
    const QMatrix4x4 tf = rgbMatrix(primaries, factor, luma);

    // linear curves, obligatory for A2B with both matrix and clut: linear-clut-curves-matrix-linear.
    cmsPipeline *a2B = cmsPipelineAlloc(nullptr, 3, 3);

    cmsToneCurve *linearCurves[3];
    linearCurves[0] = linearCurves[1] = linearCurves[2] = cmsBuildGamma(nullptr, 1.0);
    cmsStage *lCurveStage = cmsStageAllocToneCurves(nullptr, 3, linearCurves);
    cmsPipelineInsertStage(a2B, cmsAT_END, lCurveStage);

    cmsStage *clutStage = cmsStageAllocCLut16bit(nullptr, lutSize, 3, 3, nullptr);
    perceptualDummyHelper helper;
    helper.toXYZ = true;
    helper.luma = luma;
    cmsStageSampleCLut16bit(clutStage, &samplePQDummyClut, &helper, 0);
    cmsPipelineInsertStage(a2B, cmsAT_END, clutStage);

    // close to "PQ" curve
    cmsToneCurve *curve[3];
    curve[0] = curve[1] = curve[2] = cmsBuildGamma(nullptr, 2.2);
    cmsStage *curveStage = cmsStageAllocToneCurves(nullptr, 3, curve);
    cmsPipelineInsertStage(a2B, cmsAT_END, curveStage);

    // matrix

    double m[] = {
        tf(0, 0), tf(0, 1), tf(0, 2),
        tf(1, 0), tf(1, 1), tf(1, 2),
        tf(2, 0), tf(2, 1), tf(2, 2)
    };

    cmsStage* matrixStage = cmsStageAllocMatrix(nullptr, 3, 3, m, nullptr);
    cmsPipelineInsertStage(a2B, cmsAT_END, matrixStage);

    // curves, linear
    cmsPipelineInsertStage(a2B, cmsAT_END, cmsStageDup(lCurveStage));

    bool success = cmsWriteTag(iccProfile, tag, a2B);
    cmsFreeToneCurve(curve[0]);
    cmsFreeToneCurve(linearCurves[0]);
    cmsPipelineFree(a2B);
    return success;
}

bool LcmsPredefinedPipelineFunctions::setPerceptualQuantizerBToADummyPipeline(cmsHPROFILE iccProfile, cmsTagSignature tag, ColorPrimaries primaries)
{
    // From XYZ to profile space.
    const double factor = 1.0;
    const int lutSize = 4;
    cmsPipeline *b2A = cmsPipelineAlloc(nullptr, 3, 3);

           // curves, linear
    cmsToneCurve *linearCurves[3];
    linearCurves[0] = linearCurves[1] = linearCurves[2] = cmsBuildGamma(nullptr, 1.0);
    cmsStage *lCurveStage = cmsStageAllocToneCurves(nullptr, 3, linearCurves);
    cmsPipelineInsertStage(b2A, cmsAT_END, lCurveStage);

    // matrix
    double luma[3];
    const QMatrix4x4 tf = rgbMatrix(primaries, factor, luma).inverted();
    double m[] = {
        tf(0, 0), tf(0, 1), tf(0, 2),
        tf(1, 0), tf(1, 1), tf(1, 2),
        tf(2, 0), tf(2, 1), tf(2, 2)
    };

    cmsStage* matrixStage = cmsStageAllocMatrix(nullptr, 3, 3, m, nullptr);
    cmsPipelineInsertStage(b2A, cmsAT_END, matrixStage);
    // pq curve
    cmsToneCurve *curve[3];
    curve[0] = curve[1] = curve[2] = cmsBuildGamma(nullptr, 1/2.2);
    cmsStage *curveStage = cmsStageAllocToneCurves(nullptr, 3, curve);
    cmsPipelineInsertStage(b2A, cmsAT_END, curveStage);

    cmsStage *clutStage = cmsStageAllocCLut16bit(nullptr, lutSize, 3, 3, nullptr);
    perceptualDummyHelper helper;
    helper.toXYZ = false;
    helper.luma = luma;
    cmsStageSampleCLut16bit(clutStage, &samplePQDummyClut, &helper, 0);
    cmsPipelineInsertStage(b2A, cmsAT_END, clutStage);

    cmsPipelineInsertStage(b2A, cmsAT_END, cmsStageDup(lCurveStage));

    bool success = cmsWriteTag(iccProfile, tag, b2A);
    cmsFreeToneCurve(curve[0]);
    cmsFreeToneCurve(linearCurves[0]);
    cmsPipelineFree(b2A);
    return success;
}
