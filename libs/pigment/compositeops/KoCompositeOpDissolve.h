/*
 * SPDX-FileCopyrightText: 2011 Silvio Heinrich <plassy@web.de>
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
*/

#ifndef _KOCOMPOSITEOP_DISSOLVE_H_
#define _KOCOMPOSITEOP_DISSOLVE_H_

#include <KoColorSpaceMaths.h>
#include <KoCompositeOp.h>

#include <QRandomGenerator>

template<class Traits>
class KoCompositeOpDissolve: public KoCompositeOp
{
    typedef typename Traits::channels_type channels_type;

    static const qint32 channels_nb = Traits::channels_nb;
    static const qint32 alpha_pos   = Traits::alpha_pos;

    inline static quint8 getRandomValue(quint32 i)
    {
        static const quint8 randomValues[256] =
        {
            0x50, 0xAD, 0x7D, 0xA9, 0x10, 0x75, 0xCA, 0x57, 0xE2, 0x06, 0x77, 0x39, 0xD9, 0xFA, 0x5C, 0x24,
            0xEB, 0x1A, 0x6F, 0x15, 0xE7, 0x8B, 0x11, 0x71, 0xF0, 0xB9, 0x44, 0x8A, 0x27, 0x5E, 0xA1, 0x6A,
            0x47, 0x94, 0x03, 0xD5, 0xB7, 0x56, 0xEF, 0x45, 0xED, 0xBE, 0xE8, 0xB2, 0x4C, 0x0D, 0x65, 0x9E,
            0x55, 0xD7, 0x30, 0x0F, 0x52, 0xA6, 0x4D, 0x86, 0xAF, 0x66, 0x33, 0x6B, 0x3E, 0x89, 0xBD, 0xFB,
            0x00, 0xC4, 0x36, 0xFC, 0x8D, 0x4E, 0x19, 0x3F, 0x91, 0xC1, 0x40, 0x14, 0x67, 0x80, 0x17, 0x3A,
            0xF2, 0xB4, 0xD1, 0xFF, 0x35, 0xA7, 0xF7, 0x1C, 0x84, 0x2A, 0xBF, 0x46, 0xC6, 0x2B, 0x98, 0x41,
            0xF4, 0xB8, 0xA0, 0x78, 0x5A, 0xBC, 0x3B, 0x62, 0xB6, 0x7A, 0x2E, 0x07, 0x8F, 0x4A, 0xAB, 0x2F,
            0x79, 0x54, 0x81, 0x69, 0x18, 0x4F, 0xA5, 0x21, 0xD3, 0x26, 0x7C, 0x9F, 0xCF, 0xB0, 0x34, 0xCC,
            0x8C, 0xAA, 0xDB, 0x32, 0xE5, 0x1F, 0x7B, 0x37, 0x64, 0x0A, 0xF9, 0x63, 0xF5, 0x38, 0x13, 0xA2,
            0x12, 0xF6, 0xC9, 0x5D, 0xDF, 0xC7, 0x97, 0xC0, 0x51, 0xE1, 0x9A, 0x58, 0x76, 0xC3, 0x83, 0xC2,
            0x04, 0x22, 0x60, 0x9D, 0xF1, 0x5F, 0xEC, 0x6D, 0x4B, 0xE0, 0x6C, 0xD8, 0xAC, 0x25, 0xA8, 0x1E,
            0x96, 0x7E, 0x49, 0x61, 0xCD, 0x0E, 0xE3, 0xC5, 0x7F, 0x5B, 0x05, 0x6E, 0xBA, 0x0C, 0x8E, 0xF8,
            0x82, 0xDA, 0x72, 0x01, 0x23, 0x9B, 0xD2, 0x99, 0xE9, 0xC8, 0xB1, 0x28, 0xD4, 0xAE, 0x48, 0xFD,
            0x95, 0x2C, 0xE4, 0x93, 0x09, 0x3D, 0x70, 0x85, 0x43, 0x20, 0xBB, 0xDE, 0x90, 0xB3, 0x3C, 0xDD,
            0xA3, 0x73, 0x9C, 0x16, 0xDC, 0x42, 0xEA, 0x74, 0x92, 0xE6, 0xCB, 0x53, 0x08, 0xEE, 0x59, 0x02,
            0xF3, 0x29, 0xFE, 0xA4, 0x1B, 0xD6, 0x87, 0xB5, 0xCE, 0x1D, 0x68, 0x88, 0x31, 0x0B, 0x2D, 0xD0
        };

        return randomValues[i];
    }

public:
    KoCompositeOpDissolve(const KoColorSpace* cs, const QString& category)
        : KoCompositeOp(cs, COMPOSITE_DISSOLVE, category) { }

    using KoCompositeOp::composite;

    void composite(const KoCompositeOp::ParameterInfo& params) const override {

        const QBitArray& flags       = params.channelFlags.isEmpty() ? QBitArray(channels_nb,true) : params.channelFlags;
        bool             alphaLocked = (alpha_pos != -1) && !flags.testBit(alpha_pos);

        using namespace Arithmetic;

//         quint32       ctr       = quint32(reinterpret_cast<quint64>(dstRowStart) % 256);
        qint32        srcInc    = (params.srcRowStride == 0) ? 0 : channels_nb;
        bool          useMask   = params.maskRowStart != 0;
        channels_type unitValue = KoColorSpaceMathsTraits<channels_type>::unitValue;
        channels_type opacity   = KoColorSpaceMaths<float,channels_type>::scaleToA(params.opacity);

        const quint8 *srcRowStart = params.srcRowStart;
        quint8 *dstRowStart = params.dstRowStart;
        const quint8 *maskRowStart = params.maskRowStart;

        qint32 rows = params.rows;

        for(; rows>0; --rows) {
            const channels_type* src  = reinterpret_cast<const channels_type*>(srcRowStart);
            channels_type*       dst  = reinterpret_cast<channels_type*>(dstRowStart);
            const quint8*        mask = maskRowStart;

            for(qint32 c=params.cols; c>0; --c) {
                channels_type srcAlpha = (alpha_pos == -1) ? unitValue : src[alpha_pos];
                channels_type dstAlpha = (alpha_pos == -1) ? unitValue : dst[alpha_pos];
                channels_type blend    = useMask ? mul(opacity, scale<channels_type>(*mask), srcAlpha) : mul(opacity, srcAlpha);

                if (QRandomGenerator::global()->bounded(256) <= scale<quint8>(blend)
                    && blend != KoColorSpaceMathsTraits<channels_type>::zeroValue) {
                    for(qint32 i=0; i <channels_nb; i++) {
                        if(i != alpha_pos && flags.testBit(i))
                            dst[i] = src[i];
                    }

                    if(alpha_pos != -1)
                        dst[alpha_pos] = alphaLocked ? dstAlpha : unitValue;
                }

                src += srcInc;
                dst += channels_nb;
//                 ctr  = (ctr + 1) % 256;
                if (mask) {
                    ++mask;
                }
            }

            srcRowStart  += params.srcRowStride;
            dstRowStart  += params.dstRowStride;
            if (maskRowStart) {
                maskRowStart += params.maskRowStride;
            }
        }
    }
};

#endif // _KOCOMPOSITEOP_DISSOLVE_H_
