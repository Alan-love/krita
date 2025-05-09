/*
 *  SPDX-FileCopyrightText: 2017 Dmitry Kazakov <dimula73@gmail.com>
 *  SPDX-FileCopyrightText: 2022 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef KO_SVG_TEXT_SHAPE_LAYOUT_FUNC_H
#define KO_SVG_TEXT_SHAPE_LAYOUT_FUNC_H

#include "KoSvgTextShape_p.h"

namespace KoSvgTextShapeLayoutFunc
{

void calculateLineHeight(CharacterResult cr, double &ascent, double &descent, bool isHorizontal, bool compare = false);

void addWordToLine(QVector<CharacterResult> &result,
                   QPointF &currentPos,
                   QVector<int> &wordIndices,
                   LineBox &currentLine,
                   bool ltr,
                   bool isHorizontal);

void finalizeLine(QVector<CharacterResult> &result,
                  QPointF &currentPos,
                  LineBox &currentLine,
                  QPointF &lineOffset,
                  const KoSvgText::TextAnchor anchor,
                  const KoSvgText::WritingMode writingMode,
                  const bool ltr,
                  const bool inlineSize,
                  const bool textInShape,
                  const KoSvgText::ResolutionHandler &resHandler);

QVector<LineBox> breakLines(const KoSvgTextProperties &properties,
                            const QMap<int, int> &logicalToVisual,
                            QVector<CharacterResult> &result,
                            QPointF startPos, const KoSvgText::ResolutionHandler &resHandler);

QList<QPainterPath>
getShapes(QList<KoShape *> shapesInside, QList<KoShape *> shapesSubtract, const KoSvgTextProperties &properties);

QVector<LineBox> flowTextInShapes(const KoSvgTextProperties &properties,
                                  const QMap<int, int> &logicalToVisual,
                                  QVector<CharacterResult> &result,
                                  QList<QPainterPath> shapes, QPointF &startPos, const KoSvgText::ResolutionHandler &resHandler);

} // namespace KoSvgTextShapeLayoutFunc

#endif // KO_SVG_TEXT_SHAPE_LAYOUT_FUNC_H
