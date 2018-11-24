/*
 *      This is an interactive demo application for the algorithm proposed in:
 *
 *      Weighted Linde-Buzo Gray Stippling
 *      Oliver Deussen, Marc Spicker, Qian Zheng
 *
 *      In: ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia 2017)
 *      https://doi.org/10.1145/3130800.3130819
 *
 *     Copyright 2017 Marc Spicker (marc.spicker@googlemail.com)
 */

#include <QApplication>
#include <QSvgGenerator>
#include <cassert>

#include "lbgstippling.h"
#include "mainwindow.h"

QImage foveaSampling() {
    auto ellipticalGauss2DAppox = [](float x, float y, //
                                     float sigmaX, float sigmaY) -> float {
        return qExp(-((x * x) / (2.0 * sigmaX * sigmaX) + (y * y) / (2.0 * sigmaY * sigmaY)));
    };

    const QSize screenSizePx(1920, 1080);
    const QSizeF screenSizeCm(60, 33.5);
    const qreal viewDistanceCm = 80;
    const qreal foveaAlpha = 4.0 / 180.0 * M_PI;
    const qreal foveaCm = viewDistanceCm * qSin(foveaAlpha);
    const QSizeF foveaPx(screenSizePx.width() / screenSizeCm.width() * foveaCm,
                         screenSizePx.height() / screenSizeCm.height() * foveaCm);

    QImage gaussian(screenSizePx.width() * 3, screenSizePx.height() * 3, QImage::Format_Grayscale8);
    //QImage gaussian(screenSizePx.width() , screenSizePx.height() , QImage::Format_Grayscale8);
    for (int y = 0; y < gaussian.height(); ++y) {
        uchar* line = gaussian.scanLine(y);
        for (int x = 0; x < gaussian.width(); ++x) {
            float g = ellipticalGauss2DAppox(x - gaussian.width() / 2, y - gaussian.height() / 2, //
                                             foveaPx.width(), foveaPx.height());
            line[x] = qMin(static_cast<int>((1.0 - g) * 255.0), 254);
        }
    }

    return gaussian;
}
int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    app.setApplicationName("Weighted Linde-Buzo-Gray Stippling");

    QCommandLineParser parser;
    parser.setApplicationDescription("This is a helping text.");
    parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);
    parser.addHelpOption();

    const char* c = "main";
    parser.addOptions({{"input", QCoreApplication::translate(c, "Input image."),
                        QCoreApplication::translate(c, "image")},
                       {"ip", QCoreApplication::translate(c, "Initial number of points."),
                        QCoreApplication::translate(c, "points")},
                       {"ips", QCoreApplication::translate(c, "Initial point size."),
                        QCoreApplication::translate(c, "size")},
                       {"cps", QCoreApplication::translate(c, "Use constant point size.")},
                       {"psmin", QCoreApplication::translate(c, "Minimum point size."),
                        QCoreApplication::translate(c, "size")},
                       {"psmax", QCoreApplication::translate(c, "Maximum point size."),
                        QCoreApplication::translate(c, "size")},
                       {"psmapping", QCoreApplication::translate(c, "Point size mapping."),
                        QCoreApplication::translate(c, "mapping")},
                       {"ss", QCoreApplication::translate(c, "Super sampling factor."),
                        QCoreApplication::translate(c, "factor")},
                       {"iter", QCoreApplication::translate(c, "Maximum number iterations."),
                        QCoreApplication::translate(c, "number")},
                       {"hyst", QCoreApplication::translate(c, "Algorithm hysteresis."),
                        QCoreApplication::translate(c, "value")},
                       {"hystd", QCoreApplication::translate(c, "Hysteresis delta per iteration."),
                        QCoreApplication::translate(c, "value")},
                       {"output", QCoreApplication::translate(c, "Output file."),
                        QCoreApplication::translate(c, "file")},
                       {"params", QCoreApplication::translate(c, "JSON parameter file."),
                        QCoreApplication::translate(c, "params")}});

    parser.process(app);

    QImage density = foveaSampling(); // QImage(":/input/input1.jpg");
    LBGStippling::Params params;

    if (parser.isSet("input")) {
        density = QImage(parser.value("input"));
    }
    if (parser.isSet("ip")) {
        params.initialPoints = parser.value("ip").toInt();
    }
    if (parser.isSet("ips")) {
        params.initialPointSize = parser.value("ips").toFloat();
    }
    params.adaptivePointSize = !parser.isSet("cps");

    if (parser.isSet("psmin")) {
        params.pointSizeMin = parser.value("psmin").toFloat();
    }
    if (parser.isSet("psmax")) {
        params.pointSizeMax = parser.value("psmax").toFloat();
    }
    if (parser.isSet("psmapping")) {
        int index = parser.value("psmapping").toInt();
        using pmf = LBGStippling::PointMappingFunction;
        switch (index) {
        case 0:
            params.mapping = pmf::LINEAR;
            break;
        case 1:
            params.mapping = pmf::SQUAREROOT;
            break;
        case 2:
            params.mapping = pmf::EXPONENTIAL;
            break;
        case 3:
            params.mapping = pmf::SQUARE;
        }
    }
    if (parser.isSet("ss")) {
        params.superSamplingFactor = parser.value("ss").toInt();
    }
    if (parser.isSet("iter")) {
        params.maxIterations = parser.value("iter").toInt();
    }
    if (parser.isSet("hyst")) {
        params.hysteresis = parser.value("hyst").toFloat();
    }
    if (parser.isSet("hystd")) {
        params.hysteresisDelta = parser.value("hystd").toFloat();
    }
    if (parser.isSet("params")) {
        params.loadParametersJSON(parser.value("params"));
    }
    if (!parser.isSet("output")) {
        MainWindow* window = new MainWindow(density, params);
        window->show();
    } else {
        //QString outputPath = parser.value("output");
        //QStringList outputList = outputPath.split(",");
       // QStringList inputList = parser.value("input").split(",");

      //  assert(inputList.size() > outputList.size());

        LBGStippling stippling = LBGStippling();

            params.mapping =  LBGStippling::PointMappingFunction::LINEAR;

//        if (outputList.size() == inputList.size()) {
//            // for each input one output

//            for (int i = 0; i < inputList.size(); ++i) {
//                const QString& in = inputList.at(i);
//                const QString& out = outputList.at(i);
//                //density = QImage(in);
                auto result = stippling.stipple(density, params);
                std::vector<Stipple> stipples = result.stipples;
                auto map = result.indexMap;

                QImage stippleMap(stipples.size(), 2, QImage::Format_RGB32);
                QRgb* stippleLine0 = reinterpret_cast<QRgb*>(stippleMap.scanLine(0));
                QRgb* stippleLine1 = reinterpret_cast<QRgb*>(stippleMap.scanLine(1));
//                QRgb* stippleLine2 = reinterpret_cast<QRgb*>(stippleMap.scanLine(2));
//                QRgb* stippleLine3 = reinterpret_cast<QRgb*>(stippleMap.scanLine(3));
//                QRgb* stippleLine4 = reinterpret_cast<QRgb*>(stippleMap.scanLine(4));
//                QRgb* stippleLine5 = reinterpret_cast<QRgb*>(stippleMap.scanLine(5));
//                QRgb* stippleLine6 = reinterpret_cast<QRgb*>(stippleMap.scanLine(6));
//                QRgb* stippleLine7 = reinterpret_cast<QRgb*>(stippleMap.scanLine(7));
                for (int i = 0; i < stipples.size(); ++i) {
                    // Pos
                    stippleLine0[i] = static_cast<uint>(stipples[i].pos.x() * density.width() );
                    stippleLine1[i] = static_cast<uint>(stipples[i].pos.y() * density.height() );

//                    // Inner indices
//                    stippleLine2[i] = i;
//                    stippleLine3[i] = i;
//                    stippleLine4[i] = i;

//                    // Adjacent indices
//                    stippleLine5[i] = i;
//                    stippleLine6[i] = i;
//                    stippleLine7[i] = i;
                }
                stippleMap.save("stippleMap.png");

                QImage indexMap(map.width, map.height, QImage::Format_RGB32);
                for (size_t y = 0; y < map.height; ++y) {
                    QRgb* indexMapLine = reinterpret_cast<QRgb*>(indexMap.scanLine(y));
                    for (size_t x = 0; x < map.width; ++x) {
                        uint index = map.get(x, y);
                        indexMapLine[x] = index;
                    }
                }
                indexMap.save("indexMap.png");

//                QSvgGenerator generator;
//                generator.setFileName(out);
//                generator.setSize(density.size());
//                generator.setViewBox(QRectF(0, 0, density.width(), density.height()));
//                generator.setTitle("Stippling Result");
//                generator.setDescription("SVG File created by Weighted Linde-Buzo-Gray Stippling");

//                QPainter painter;
//                painter.begin(&generator);
//                painter.setPen(Qt::NoPen);

//                for (const auto& s : stipples) {
//                    auto x = s.pos.x() * density.width();
//                    auto y = s.pos.y() * density.height();
//                    painter.setBrush(s.color);
//                    painter.drawEllipse(QPointF(x, y), s.size / 2.0, s.size / 2.0);
//                }

//                painter.end();
//            }
//        } else if (outputList.size() == 1) {
//            // merge all inputs into one output

//            QSvgGenerator generator;
//            generator.setFileName(outputList.at(0));
//            generator.setSize(density.size());
//            generator.setViewBox(QRectF(0, 0, density.width(), density.height()));
//            generator.setTitle("Stippling Result");
//            generator.setDescription("SVG File created by Weighted Linde-Buzo-Gray Stippling");

//            QPainter painter;
//            painter.begin(&generator);
//            painter.setPen(Qt::NoPen);

//            for (const QString& in : inputList) {
//                density = QImage(in);
//                std::vector<Stipple> stipples = stippling.stipple(density, params).stipples;

//                for (const auto& s : stipples) {
//                    auto x = s.pos.x() * density.width();
//                    auto y = s.pos.y() * density.height();
//                    painter.setBrush(s.color);
//                    painter.drawEllipse(QPointF(x, y), s.size / 2.0, s.size / 2.0);
//                }
//            }
//            painter.end();
//        }
        exit(0);
    }

    return app.exec();
}
