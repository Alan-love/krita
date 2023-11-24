/*
 * SPDX-FileCopyrightText: 2023 Alvin Wong <alvin@alvinhc.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SVG_CREATE_TEXT_STRATEGY_H
#define SVG_CREATE_TEXT_STRATEGY_H

#include <KoInteractionStrategy.h>

#include <QPointF>

class SvgTextTool;

class KoSvgTextShape;

class SvgCreateTextStrategy : public KoInteractionStrategy
{
public:
    SvgCreateTextStrategy(SvgTextTool *tool, const QPointF &clicked, double pressure = 0.0);
    ~SvgCreateTextStrategy() override = default;

    void paint(QPainter &painter, const KoViewConverter &converter) override;
    void handleMouseMove(const QPointF &mouseLocation, Qt::KeyboardModifiers modifiers) override;
    KUndo2Command *createCommand() override;
    void cancelInteraction() override;
    void finishInteraction(Qt::KeyboardModifiers modifiers) override;

private:
    QPointF m_dragStart;
    QPointF m_dragEnd;
    double m_pressure;
    Qt::KeyboardModifiers m_modifiers;
};

#endif /* SVG_CREATE_TEXT_STRATEGY_H */
