#-------------------------------------------------
#
# Project created by QtCreator 2016-04-01T09:59:44
#
#-------------------------------------------------

QT       -= core
QT       -= gui

TARGET = SuGeT
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    trees/genespeciestreeutil.cpp \
    trees/newicklex.cpp \
    trees/node.cpp \
    trees/treeinfo.cpp \
    trees/treeiterator.cpp \
    supergenetreemaker.cpp \
    genesubtreecorrector.cpp \
    trees/polysolver.cpp


QMAKE_CXXFLAGS += -std=c++0x

HEADERS += \
    trees/genespeciestreeutil.h \
    trees/newicklex.h \
    trees/node.h \
    trees/treeinfo.h \
    trees/treeiterator.h \
    supergenetreemaker.h \
    div/util.h \
    genesubtreecorrector.h \
    div/tinydir.h \
    trees/polysolver.h


DEFINES -= UNICODE
DEFINES += _MBCS
