/*
 *  SPDX-FileCopyrightText: 2026 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "LcmsPredefinedPipelineFunctions.h"
#include <KoColorTransferFunctions.h>
#include <KoColorProfile.h>
#include <QTransform>
#include <QGenericMatrix>

/**
 * @brief bradfordMatrix
 * From http://brucelindbloom.com/index.html?Eqn_ChromAdapt.html
 * @param src -- source white point XYZ
 * @param dst -- destionation white point XYZ
 * @return chromatic adaptation matrix.
 */
static QTransform bradfordMatrix(double *src, double *dst) {
    double bradford[] = {
        0.8951000, 0.2664000, -0.1614000,
        -0.7502000, 1.7135000, 0.0367000,
        0.0389000, -0.0685000, 1.0296000
    };

    QGenericMatrix<1, 3, double> srcLMS = QGenericMatrix<3, 3, double>(bradford) * QGenericMatrix<1, 3, double>(src);
    QGenericMatrix<1, 3, double> dstLMS = QGenericMatrix<3, 3, double>(bradford) * QGenericMatrix<1, 3, double>(dst);

    QTransform adaptation (
        dstLMS(0, 0)/srcLMS(0, 0), 0, 0,
        0, dstLMS(1, 0)/srcLMS(1, 0), 0,
        0, 0, dstLMS(2, 0)/srcLMS(2, 0)
    );

    QTransform bradfordInv (
        0.9869929, -0.1470543, 0.1599627,
        0.4323053, 0.5183603, 0.0492912,
        -0.0085287, 0.0400428, 0.9684867
        );

    return bradfordInv * adaptation * bradfordInv.inverted();
}
/**
 * @brief rgbMatrix
 * Generates a proper rgb matrix from a primaries type.
 * @param primaries -- default primaries to use.
 * @param factor -- scaling factor.
 * @return transfrom representing a XYZ to RGB conversion.
 */
static QTransform rgbMatrix(ColorPrimaries primaries, double factor = 1.0) {
    const double mul = 1.0/factor;
    // Calculation: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    // Calculates XYZ with Y = 1.0
    QVector<double> col;
    KoColorProfile::colorantsForType(primaries, col);
    QTransform initialXYZMatrix(
        col[2]/col[3], col[4]/col[5], col[6]/col[7],
        1.0, 1.0, 1.0,
        (1.0-col[2]-col[3])/col[3], (1.0-col[4]-col[5])/col[5], (1.0-col[6]-col[7])/col[7]
        );
    QTransform invTf = initialXYZMatrix.inverted();

    // Calculate scaling factor with white point.
    double wp [] = {
        col[0]/col[1], 1.0, (1.0-col[0]-col[1])/col[1]
    };
    double inv[] = {
        invTf.m11(), invTf.m12(), invTf.m13(),
        invTf.m21(), invTf.m22(), invTf.m23(),
        invTf.m31(), invTf.m32(), invTf.m33()
    };

    QGenericMatrix<1, 3, double> S = QGenericMatrix<3, 3, double>(inv) * QGenericMatrix<1, 3, double>(wp);

    // Apply scaling factor to normalize the Y and get the matrix.
    QTransform finalXYZMatrix(
        initialXYZMatrix.m11() * S(0, 0), initialXYZMatrix.m12() * S(1, 0), initialXYZMatrix.m13() * S(2, 0),
        initialXYZMatrix.m21() * S(0, 0), initialXYZMatrix.m22() * S(1, 0), initialXYZMatrix.m23() * S(2, 0),
        initialXYZMatrix.m31() * S(0, 0), initialXYZMatrix.m32() * S(1, 0), initialXYZMatrix.m33() * S(2, 0)
        );
    // Apply bradford white point adaptation.
    double wpD50 [] = {
        0.9642, 1.0, 0.8249
    };
    finalXYZMatrix = bradfordMatrix(wp, wpD50) * finalXYZMatrix;

    return QTransform(
        finalXYZMatrix.m11() * mul, finalXYZMatrix.m12() * mul, finalXYZMatrix.m13() * mul,
        finalXYZMatrix.m21() * mul, finalXYZMatrix.m22() * mul, finalXYZMatrix.m23() * mul,
        finalXYZMatrix.m31() * mul, finalXYZMatrix.m32() * mul, finalXYZMatrix.m33() * mul
        );
}

struct perceptualDummyHelper {
    bool toXYZ{false};
    int channels {3};
    };

/**
 * @brief samplePQDummyClut
 * In this sampler, we convert the PQ data to linear and then to HLG. The benefit of HLG
 * is that it sort of allows for tonemapping.
 */
cmsInt32Number samplePQDummyClut(const cmsUInt16Number In[], cmsUInt16Number Out[], void *Cargo) {
    struct perceptualDummyHelper *helper = (struct perceptualDummyHelper *) Cargo;
    const float scale = 3.7743;
    for (int i = 0; i < helper->channels; i++) {
        const float val = float(In[i])/65535.0;

        cmsUInt16Number finalResult = In[i];
        if (!helper->toXYZ) {
            const float lin = removeSmpte2048Curve(val);
            finalResult = qBound(0, qRound( applyHLGCurve(lin / scale) * 65535 ), 65535);
        } else {
            const float lin = removeHLGCurve(val) * scale;
            finalResult = qBound(0, qRound( applySmpte2048Curve(lin) * 65535 ), 65535);
        }
        Out[i] = finalResult;
    }
    return true;
};

bool LcmsPredefinedPipelineFunctions::setPerceptualQuantizerAToBDummyPipeline(cmsHPROFILE iccProfile, cmsTagSignature tag, ColorPrimaries primaries)
{
    // From profile space to XYZ.
    const double factor = 2.0;



    // linear curves, obligatory for A2B with both matrix and clut: linear-clut-curves-matrix-linear.
    cmsPipeline *a2B = cmsPipelineAlloc(nullptr, 3, 3);

    cmsToneCurve *linearCurves[3];
    linearCurves[0] = linearCurves[1] = linearCurves[2] = cmsBuildGamma(nullptr, 1.0);
    cmsStage *lCurveStage = cmsStageAllocToneCurves(nullptr, 3, linearCurves);
    cmsPipelineInsertStage(a2B, cmsAT_END, lCurveStage);

    cmsStage *clutStage = cmsStageAllocCLut16bit(nullptr, 8, 3, 3, nullptr);
    perceptualDummyHelper helper;
    cmsStageSampleCLut16bit(clutStage, &samplePQDummyClut, &helper, 0);
    cmsPipelineInsertStage(a2B, cmsAT_END, clutStage);

    // close to "PQ" curve
    cmsToneCurve *curve[3];
    curve[0] = curve[1] = curve[2] = cmsBuildGamma(nullptr, 2.2);
    cmsStage *curveStage = cmsStageAllocToneCurves(nullptr, 3, curve);
    cmsPipelineInsertStage(a2B, cmsAT_END, curveStage);

    // matrix

    QTransform tf = rgbMatrix(primaries, factor);
    double m[] = {
        tf.m11(), tf.m12(), tf.m13(),
        tf.m21(), tf.m22(), tf.m23(),
        tf.m31(), tf.m32(), tf.m33()
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
    cmsPipeline *b2A = cmsPipelineAlloc(nullptr, 3, 3);

           // curves, linear
    cmsToneCurve *linearCurves[3];
    linearCurves[0] = linearCurves[1] = linearCurves[2] = cmsBuildGamma(nullptr, 1.0);
    cmsStage *lCurveStage = cmsStageAllocToneCurves(nullptr, 3, linearCurves);
    cmsPipelineInsertStage(b2A, cmsAT_END, lCurveStage);

    // matrix
    QTransform tf = rgbMatrix(primaries, 2.0).inverted();
    double m[] = {
        tf.m11(), tf.m12(), tf.m13(),
        tf.m21(), tf.m22(), tf.m23(),
        tf.m31(), tf.m32(), tf.m33()
    };

    cmsStage* matrixStage = cmsStageAllocMatrix(nullptr, 3, 3, m, nullptr);
    cmsPipelineInsertStage(b2A, cmsAT_END, matrixStage);
    // pq curve
    cmsToneCurve *curve[3];
    curve[0] = curve[1] = curve[2] = cmsBuildGamma(nullptr, 1/2.2);
    cmsStage *curveStage = cmsStageAllocToneCurves(nullptr, 3, curve);
    cmsPipelineInsertStage(b2A, cmsAT_END, curveStage);

    cmsStage *clutStage = cmsStageAllocCLut16bit(nullptr, 8, 3, 3, nullptr);
    perceptualDummyHelper helper;
    helper.toXYZ = false;
    cmsStageSampleCLut16bit(clutStage, &samplePQDummyClut, &helper, 0);
    cmsPipelineInsertStage(b2A, cmsAT_END, clutStage);

    cmsPipelineInsertStage(b2A, cmsAT_END, cmsStageDup(lCurveStage));

    bool success = cmsWriteTag(iccProfile, tag, b2A);
    cmsFreeToneCurve(curve[0]);
    cmsFreeToneCurve(linearCurves[0]);
    cmsPipelineFree(b2A);
    return success;
}
