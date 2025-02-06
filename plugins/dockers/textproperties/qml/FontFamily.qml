/*
 *  SPDX-FileCopyrightText: 2024 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
import QtQuick 2.0
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.12
import org.krita.flake.text 1.0

CollapsibleGroupProperty {
    propertyName: i18nc("@label", "Font Family");
    propertyType: TextPropertyBase.Character;
    toolTip: i18nc("@info:tooltip",
                   "Font family allows selecting a list of fonts that should be used for the current text. The first font family is the primary font used, while each font family after that is used for fallback.");
    searchTerms: i18nc("comma separated search terms for the font-family property, matching is case-insensitive",
                       "font-family, typeface, font, fallback");

    property var fontFamilies: [];
    property var fontFamilyModel : [];

    onPropertiesUpdated: {
        blockSignals = true;
        fontFamilies = properties.fontFamilies;
        visible = properties.fontFamiliesState !== KoSvgTextPropertiesModel.PropertyUnset;
        blockSignals = false;
    }

    onFontFamiliesChanged: {
        if (blockSignals) {
            if (fontFamilies.length >0) {
                mainFamilyCmb.currentIndex = mainFamilyCmb.find(fontFamilies[0]);
            }
            familyListView.model = fontFamilies;
        } else {
            properties.fontFamilies = fontFamilies;
        }
        mainWindow.slotUpdateStylesModel();
    }

    titleItem: RowLayout {
        width: parent.width;
        height: childrenRect.height;
        Label {
            id: propertyTitle;
            text: propertyName;
            verticalAlignment: Text.AlignVCenter
            color: sysPalette.text;
            elide: Text.ElideRight;
            Layout.maximumWidth: contentWidth;

        }
        ComboBox {
        id: mainFamilyCmb;
        model: fontFamilyModel;
        Layout.fillWidth: true;
        onActivated: {if (fontFamilies.length >0) {
                fontFamilies[0] = currentText;
            }
        }
        textRole: "name";

        delegate: ItemDelegate {
            id: fontDelegate
            required property var model;

            highlighted: mainFamilyCmb.highlightedIndex === model.index;
            contentItem: KoShapeQtQuickLabel{
                id: fontFamilyDelegate;
                property bool colorBitmap : model.metadata["color_bitmap"];
                property bool colorCLRV0 : model.metadata["color_clrv0"];
                property bool colorCLRV1 : model.metadata["color_clrv1"];
                property bool colorSVG : model.metadata["color_svg"];
                property bool isVariable : model.metadata["is_variable"];
                property int type : model.metadata["font_type"];

                Layout.fillWidth: true;
                svgData: model.metadata["sample_svg"];
                imageScale: 3;
                imagePadding: nameLabel.height;
                foregroundColor: sysPalette.text;
                fullColor: colorBitmap || colorCLRV0 || colorCLRV1 || colorSVG;

                Label {
                    id: nameLabel;
                    text: model.name;
                    anchors.top: parent.top;
                    anchors.left: parent.left;
                }

                Row {
                    id: featureRow;
                    spacing: columnSpacing;
                    anchors.bottom: parent.bottom;
                    anchors.left: parent.left;
                    property int imgHeight: 12;
                    Image {
                        property int type: fontFamilyDelegate.type;
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        source: type === KoSvgText.BDFFontType? "qrc:///light_select-pixel.svg"
                                         : type === KoSvgText.Type1FontType? "qrc:///light_transparency-enabled.svg"
                                         : type === KoSvgText.OpenTypeFontType? "qrc:///light_select-shape.svg":"qrc:///light_system-help.svg";
                    }
                    Image {
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        visible: fontFamilyDelegate.isVariable;
                        source: "qrc:///light_zoom-horizontal.svg"
                    }

                    Image {
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        visible: fontFamilyDelegate.colorBitmap;
                        source: "qrc:///light_palette-library.svg"
                    }
                    Image {
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        visible: fontFamilyDelegate.colorCLRV0;
                        source: "qrc:///light_color-adjustment-mode-channels.svg"
                    }
                    Image {
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        visible: fontFamilyDelegate.colorCLRV1;
                        source: "qrc:///light_config-color-manage.svg"
                    }
                    Image {
                        width: parent.imgHeight;
                        height: parent.imgHeight;
                        visible: fontFamilyDelegate.colorSVG;
                        source: "qrc:///light_sort-by-hue.svg"
                    }
                }
            }
            width: fontResourceView.listWidth;
            background: Rectangle { color: highlighted? parent.palette.highlight:"transparent"; }
        }
        popup: Popup {
            y: mainFamilyCmb.height - 1;
            x: mainFamilyCmb.width - width;
            width: contentWidth;
            height: contentHeight;
            padding: 2;

            contentItem: FontResourceView {
                id: fontResourceView;
                fontModel: mainFamilyCmb.delegateModel;
                tagModel: fontTagModel;
                currentIndex: mainFamilyCmb.highlightedIndex;
            }
            palette: mainFamilyCmb.palette;
        }
        wheelEnabled: true;
    }
    }

    onEnableProperty: properties.fontFamiliesState = KoSvgTextPropertiesModel.PropertySet;

    contentItem: GridLayout {
        columns: 2
        anchors.left: parent.left
        anchors.right: parent.right
        columnSpacing: columnSpacing;

        RevertPropertyButton {
            revertState: properties.fontFamiliesState;
            onClicked: properties.fontFamiliesState = KoSvgTextPropertiesModel.PropertyUnset;
        }

        ScrollView {
            id: fullFamilyList;
            Layout.fillWidth: true;
            Layout.preferredHeight: ItemDelegate.implicitHeight * 3;
            background: Rectangle {
                color: sysPalette.alternateBase;
                border.color: sysPalette.text;
                border.width: 1;
            }
            ListView {
                id: familyListView;
                anchors.fill: parent;
                model: []
                delegate: ItemDelegate {
                    text: modelData;
                    width: parent.width;
                }
                Label {
                    text: i18n("Family list is empty");
                    wrapMode: Text.WordWrap;
                    anchors.fill: parent;
                    anchors.horizontalCenter: parent.horizontalCenter;
                    visible: parent.count === 0;
                }
            }
        }


    }
}

