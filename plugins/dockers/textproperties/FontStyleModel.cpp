/*
 *  SPDX-FileCopyrightText: 2024 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#include "FontStyleModel.h"
#include <QDebug>

const QString WEIGHT_TAG = "wght";
const QString WIDTH_TAG = "wdth";
const QString SLANT_TAG = "slnt";
const QString ITALIC_TAG = "ital";
const QString OPTICAL_TAG = "opsz";

struct FontStyleModel::Private {
    QList<KoSvgText::FontFamilyStyleInfo> styles;
    QList<QLocale> locales;
};

FontStyleModel::FontStyleModel(QObject *parent)
    : QAbstractItemModel(parent)
    , d(new Private)
{

}

FontStyleModel::~FontStyleModel()
{
}
// sort the font styles in any given manner.
static bool styleLowerThan(const KoSvgText::FontFamilyStyleInfo &a, const KoSvgText::FontFamilyStyleInfo &b) {

    for (auto it = a.instanceCoords.begin(); it != a.instanceCoords.end(); it++) {
        if (it.key() == ITALIC_TAG || it.key() == SLANT_TAG) continue;
        qreal other = b.instanceCoords.value(it.key());
        if (it.value() < other) return true;
    }
    // todo: support slant and ital variation axes.
    QFont::Style aStyle = a.isItalic? a.isOblique?  QFont::StyleOblique: QFont::StyleItalic: QFont::StyleNormal;
    QFont::Style bStyle = b.isItalic? a.isOblique?  QFont::StyleOblique: QFont::StyleItalic: QFont::StyleNormal;
    return aStyle < bStyle;
}

void FontStyleModel::setStylesInfo(QList<KoSvgText::FontFamilyStyleInfo> styles)
{
    beginResetModel();
    std::sort(styles.begin(), styles.end(), styleLowerThan);
    d->styles = styles;
    endResetModel();
}

void FontStyleModel::setLocales(QList<QLocale> locales)
{
    d->locales = locales;
}

qreal FontStyleModel::weightValue(int row)
{
    return d->styles.value(row).instanceCoords.value(WEIGHT_TAG, 400.0);
}

qreal FontStyleModel::widthValue(int row)
{
    return d->styles.value(row).instanceCoords.value(WIDTH_TAG, 100.0);
}

int FontStyleModel::styleModeValue(int row)
{
    KoSvgText::FontFamilyStyleInfo style = d->styles.value(row);
    QFont::Style styleType = style.isItalic? style.isOblique? QFont::StyleOblique: QFont::StyleItalic: QFont::StyleNormal;
    if (style.instanceCoords.value(ITALIC_TAG, 0) != 0) {
        styleType = QFont::StyleItalic;
    } else if (style.instanceCoords.value(SLANT_TAG, 0) != 0) {
        styleType = QFont::StyleOblique;
    }
    return styleType;
}

qreal FontStyleModel::slantValue(int row)
{
    return -(d->styles.value(row).instanceCoords.value(SLANT_TAG, 0));
}

QVariantHash FontStyleModel::axesValues(int row)
{
    return data(index(row, 0), AxisValues).toHash();
}

// TODO: consider somehow using the exact same method as in KoFFWWSConverter
QHash<int, KoSvgText::FontFamilyStyleInfo> searchAxisTag(QString tag, qreal value, QVector<qreal> values, qreal defaultVal, QHash<int, KoSvgText::FontFamilyStyleInfo> styles) {
    QHash<int, KoSvgText::FontFamilyStyleInfo> filteredStyles;
    auto valIt = std::lower_bound(values.begin(), values.end(), value);
    if (valIt == values.end()) valIt--;
    qreal selectedValue = *valIt;
    if (valIt != values.begin()) valIt--;
    if (qAbs(*valIt - value) < qAbs(selectedValue - value)) {
        selectedValue = *valIt;
    }
    for (auto it = styles.begin(); it!= styles.end(); it++) {
        if (it.value().instanceCoords.value(tag, defaultVal) == selectedValue) {
            filteredStyles.insert(it.key(), it.value());
        }
    }
    return filteredStyles;
}

int FontStyleModel::rowForStyle(qreal weight, qreal width, int styleMode)
{
    QHash<int, KoSvgText::FontFamilyStyleInfo> styles;
    QVector<qreal> weights;
    QVector<qreal> widths;
    for (int i = 0; i< d->styles.size(); i++) {
        styles.insert(i, d->styles.at(i));
        weights.append(d->styles.value(i).instanceCoords.value(WEIGHT_TAG, 400.0));
        widths.append(d->styles.value(i).instanceCoords.value(WIDTH_TAG, 100.0));
    }
    std::sort(weights.begin(), weights.end());
    std::sort(widths.begin(), widths.end());

    if (styles.size() > 1) {
        styles = searchAxisTag(WEIGHT_TAG, weight, weights, 400.0, styles);
    }
    if (styles.size() > 1) {
        styles = searchAxisTag(WIDTH_TAG, width, widths, 100.0, styles);
    }
    if (styles.size() > 1) {
        QHash<int, KoSvgText::FontFamilyStyleInfo> filteredStyles;
        for (auto it = styles.begin(); it!= styles.end(); it++) {
            int fontStyle = it.value().isItalic? it.value().isOblique? QFont::StyleOblique: QFont::StyleItalic: QFont::StyleNormal;
            if (fontStyle == styleMode) {
                filteredStyles.insert(it.key(), it.value());
            }
        }
        if (!filteredStyles.isEmpty()) styles = filteredStyles;
    }
    return styles.isEmpty()? 0: styles.keys().first();
}

QModelIndex FontStyleModel::index(int row, int column, const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    if (column != 0) return QModelIndex();
    if (row >= 0 && row < d->styles.size()) return createIndex(row, column, &row);
    return QModelIndex();
}

QModelIndex FontStyleModel::parent(const QModelIndex &child) const
{
    Q_UNUSED(child)
     return QModelIndex();
}

int FontStyleModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return d->styles.size();
}

int FontStyleModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return 1;
}

QVariant FontStyleModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();

    KoSvgText::FontFamilyStyleInfo style = d->styles.at(index.row());
    if (role == Qt::DisplayRole) {
        if (style.localizedLabels.isEmpty()) return QString();
        QString label = style.localizedLabels.value(QLocale(QLocale::English), style.localizedLabels.values().first());
        Q_FOREACH(const QLocale &locale, d->locales) {
            if (style.localizedLabels.keys().contains(locale)) {
                label = style.localizedLabels.value(locale, label);
                break;
            }
        }
        return label;
    } else if (role == Weight){
        return style.instanceCoords.value(WEIGHT_TAG, 400.0);
    } else if (role == Width){
        return style.instanceCoords.value(WIDTH_TAG, 100.0);
    } else if (role == StyleMode){
        QFont::Style styleType = style.isItalic? style.isOblique? QFont::StyleOblique: QFont::StyleItalic: QFont::StyleNormal;
        if (style.instanceCoords.value(ITALIC_TAG, 0) != 0) {
            styleType = QFont::StyleItalic;
        } else if (style.instanceCoords.value(SLANT_TAG, 0) != 0) {
            styleType = QFont::StyleOblique;
        }
        return styleType;
    } else if (role == Slant){
        return -(style.instanceCoords.value(SLANT_TAG, 0));
    } else if (role == AxisValues) {
        QVariantHash vals;
        for (auto it = style.instanceCoords.begin(); it != style.instanceCoords.end(); it++) {
            if (!(it.key() == WEIGHT_TAG
                  || it.key() == WIDTH_TAG
                  || it.key() == SLANT_TAG
                  || it.key() == ITALIC_TAG)) {
                vals.insert(it.key(), QVariant::fromValue(it.value()));
            }
        }
        return vals;
    }
    return QVariant();
}

QHash<int, QByteArray> FontStyleModel::roleNames() const
{
    QHash<int, QByteArray> roles = QAbstractItemModel::roleNames();
    roles[Weight] = "weight";
    roles[Width] = "width";
    roles[StyleMode] = "stylemode";
    roles[Slant] = "slant";
    roles[AxisValues] = "axisvalues";
    return roles;
}
