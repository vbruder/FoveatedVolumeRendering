INCLUDEPATH += $$PWD

qtHaveModule(opengl)|qtConfig(opengles2)  {
    DEFINES += QT_OPENGL_SUPPORT
    QT += opengl widgets
}

SOURCES += \
    $$PWD/hoverpoints.cpp

HEADERS += \
    $$PWD/hoverpoints.h

#RESOURCES += $$PWD/shared.qrc

