/*
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_rectangle_constraint_widget.h"
#include "kis_tool_rectangle_base.h"

#include <kis_icon.h>
#include "kis_aspect_ratio_locker.h"
#include "kis_signals_blocker.h"
#include <KConfigGroup>
#include <KSharedConfig>

KisRectangleConstraintWidget::KisRectangleConstraintWidget(QWidget *parent, KisToolRectangleBase *tool, bool showRoundCornersGUI)
    : QWidget(parent)
{
    m_tool = tool;

    setupUi(this);

    connect(lockWidthButton, SIGNAL(toggled(bool)), this, SLOT(inputsChanged()));
    connect(lockHeightButton, SIGNAL(toggled(bool)), this, SLOT(inputsChanged()));
    connect(lockRatioButton, SIGNAL(toggled(bool)), this, SLOT(inputsChanged()));

    QIcon lockedIcon = KisIconUtils::loadIcon("locked");
    QIcon unlockedIcon = KisIconUtils::loadIcon("unlocked");
    lockWidthButton->setIcon(lockWidthButton->isChecked() ? lockedIcon : unlockedIcon);
    lockHeightButton->setIcon(lockHeightButton->isChecked() ? lockedIcon : unlockedIcon);
    lockRatioButton->setIcon(lockRatioButton->isChecked() ? lockedIcon : unlockedIcon);

    connect(intWidth,  SIGNAL(valueChanged(int)), this, SLOT(inputWidthChanged()));
    connect(intHeight, SIGNAL(valueChanged(int)), this, SLOT(inputHeightChanged()));
    connect(doubleRatio, SIGNAL(valueChanged(double)), this, SLOT(inputRatioChanged()));

    connect(this, SIGNAL(constraintsChanged(bool,bool,bool,float,float,float)), m_tool, SLOT(constraintsChanged(bool,bool,bool,float,float,float)));
    connect(m_tool, SIGNAL(rectangleChanged(QRectF)), this, SLOT(rectangleChanged(QRectF)));

    m_cornersAspectLocker = new KisAspectRatioLocker(this);
    m_cornersAspectLocker->connectSpinBoxes(intRoundCornersX, intRoundCornersY, cornersAspectButton);

    connect(m_cornersAspectLocker, SIGNAL(sliderValueChanged()), SLOT(slotRoundCornersChanged()));
    connect(m_cornersAspectLocker, SIGNAL(aspectButtonChanged()), SLOT(slotRoundCornersAspectLockChanged()));

    connect(m_tool, SIGNAL(sigRequestReloadConfig()), SLOT(slotReloadConfig()));
    slotReloadConfig();

    if (!showRoundCornersGUI) {
        intRoundCornersX->setVisible(false);
        intRoundCornersY->setVisible(false);
        cornersAspectButton->setVisible(false);
    }
}

void KisRectangleConstraintWidget::inputWidthChanged()
{
    lockWidthButton->setChecked(true);
    inputsChanged();
}

void KisRectangleConstraintWidget::inputHeightChanged()
{
    lockHeightButton->setChecked(true);
    inputsChanged();
}

void KisRectangleConstraintWidget::inputRatioChanged()
{
    lockRatioButton->setChecked(true);
    inputsChanged();
}

void KisRectangleConstraintWidget::inputsChanged()
{
    Q_EMIT constraintsChanged(
        lockRatioButton->isChecked(),
        lockWidthButton->isChecked(),
        lockHeightButton->isChecked(),
        doubleRatio->value(),
        intWidth->value(),
        intHeight->value()
                );
    QIcon lockedIcon = KisIconUtils::loadIcon("locked");
    QIcon unlockedIcon = KisIconUtils::loadIcon("unlocked");
    lockWidthButton->setIcon(lockWidthButton->isChecked() ? lockedIcon : unlockedIcon);
    lockHeightButton->setIcon(lockHeightButton->isChecked() ? lockedIcon : unlockedIcon);
    lockRatioButton->setIcon(lockRatioButton->isChecked() ? lockedIcon : unlockedIcon);
}

void KisRectangleConstraintWidget::slotRoundCornersChanged()
{
    m_tool->roundCornersChanged(intRoundCornersX->value(), intRoundCornersY->value());

    KConfigGroup cfg = KSharedConfig::openConfig()->group(m_tool->toolId());
    cfg.writeEntry("roundCornersX", intRoundCornersX->value());
    cfg.writeEntry("roundCornersY", intRoundCornersY->value());
}

void KisRectangleConstraintWidget::slotRoundCornersAspectLockChanged()
{
    KConfigGroup cfg = KSharedConfig::openConfig()->group(m_tool->toolId());
    cfg.writeEntry("roundCornersAspectLocked", cornersAspectButton->keepAspectRatio());
}

void KisRectangleConstraintWidget::slotReloadConfig()
{
    KConfigGroup cfg = KSharedConfig::openConfig()->group(m_tool->toolId());

    {
        KisSignalsBlocker b(intRoundCornersX, intRoundCornersY, cornersAspectButton);
        intRoundCornersX->setValue(cfg.readEntry("roundCornersX", 0));
        intRoundCornersY->setValue(cfg.readEntry("roundCornersY", 0));
        cornersAspectButton->setKeepAspectRatio(cfg.readEntry("roundCornersAspectLocked", true));
        m_cornersAspectLocker->updateAspect();
    }

    slotRoundCornersChanged();
}

void KisRectangleConstraintWidget::rectangleChanged(const QRectF &rect) 
{
    intWidth->blockSignals(true);
    intHeight->blockSignals(true);
    doubleRatio->blockSignals(true);

    if (!lockWidthButton->isChecked()) intWidth->setValue(rect.width());
    if (!lockHeightButton->isChecked()) intHeight->setValue(rect.height());

    if (!lockRatioButton->isChecked() && !(rect.width() == 0 && rect.height() == 0)) {
        doubleRatio->setValue(fabs(rect.width()) / fabs(rect.height()));
    }

    intWidth->blockSignals(false);
    intHeight->blockSignals(false);
    doubleRatio->blockSignals(false);
}
