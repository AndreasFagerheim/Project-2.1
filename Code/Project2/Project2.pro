TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        functions.cpp \
        main.cpp \
        test.cpp

INCLUDEPATH += /usr/local/Cellar/armadillo/9.700.2/include/
LIBS += -L/usr/local/Cellar/armadillo/9.700.2/lib/ -larmadillo

HEADERS += \
    functions.h \
    test.h

