// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

// $Id: lagrangecoeff9.cc,v 1.1 2000/03/11 14:14:00 heuvelin Exp $

/// \author Martin Baumann, Michael Schick

#include "polynomials/lagrangecoeff9.h"
#include "polynomials/lagrangecoeff.h"

namespace hiflow {

static double _poly_lagrange_coeff9[] = {1.000000000000000000000000000000e+00,
                                         -2.546071428571428571428571428570e+01,
                                         2.617633928571428571428571428570e+02,
                                         -1.453821428571428571428571428570e+03,
                                         4.869492187500000000000000000000e+03,
                                         -1.029598593750000000000000000000e+04,
                                         1.383960937500000000000000000000e+04,
                                         -1.146710491071428571428571428570e+04,
                                         5.338135044642857142857142857140e+03,
                                         -1.067627008928571428571428571430e+03,

                                         0.000000000000000000000000000000e-01,
                                         8.100000000000000000000000000000e+01,
                                         -1.333317857142857142857142857140e+03,
                                         9.202974107142857142857142857140e+03,
                                         -3.493276875000000000000000000000e+04,
                                         8.003394843750000000000000000000e+04,
                                         -1.136693250000000000000000000000e+05,
                                         9.798443437500000000000000000000e+04,
                                         -4.697558839285714285714285714290e+04,
                                         9.608643080357142857142857142860e+03,

                                         0.000000000000000000000000000000e-01,
                                         -1.620000000000000000000000000000e+02,
                                         3.395635714285714285714285714290e+03,
                                         -2.712530892857142857142857142860e+04,
                                         1.134551812500000000000000000000e+05,
                                         -2.783094187500000000000000000000e+05,
                                         4.155573375000000000000000000000e+05,
                                         -3.720087000000000000000000000000e+05,
                                         1.836318455357142857142857142860e+05,
                                         -3.843457232142857142857142857140e+04,

                                         0.000000000000000000000000000000e-01,
                                         2.520000000000000000000000000000e+02,
                                         -5.660100000000000000000000000000e+03,
                                         4.898407500000000000000000000000e+04,
                                         -2.194107750000000000000000000000e+05,
                                         5.688797062500000000000000000000e+05,
                                         -8.879493375000000000000000000000e+05,
                                         8.237335500000000000000000000000e+05,
                                         -4.185097875000000000000000000000e+05,
                                         8.968066875000000000000000000000e+04,

                                         0.000000000000000000000000000000e-01,
                                         -2.835000000000000000000000000000e+02,
                                         6.580237500000000000000000000000e+03,
                                         -5.940438750000000000000000000000e+04,
                                         2.784985031250000000000000000000e+05,
                                         -7.538794031250000000000000000000e+05,
                                         1.222683356250000000000000000000e+06,
                                         -1.172491706250000000000000000000e+06,
                                         6.128179031250000000000000000000e+05,
                                         -1.345210031250000000000000000000e+05,

                                         0.000000000000000000000000000000e-01,
                                         2.268000000000000000000000000000e+02,
                                         -5.366250000000000000000000000000e+03,
                                         4.970868750000000000000000000000e+04,
                                         -2.402510625000000000000000000000e+05,
                                         6.719489156250000000000000000000e+05,
                                         -1.125621562500000000000000000000e+06,
                                         1.112704593750000000000000000000e+06,
                                         -5.978711250000000000000000000000e+05,
                                         1.345210031250000000000000000000e+05,

                                         0.000000000000000000000000000000e-01,
                                         -1.260000000000000000000000000000e+02,
                                         3.019050000000000000000000000000e+03,
                                         -2.845361250000000000000000000000e+04,
                                         1.405010812500000000000000000000e+05,
                                         -4.028043937500000000000000000000e+05,
                                         6.930876375000000000000000000000e+05,
                                         -7.041593250000000000000000000000e+05,
                                         3.886162312500000000000000000000e+05,
                                         -8.968066875000000000000000000000e+04,

                                         0.000000000000000000000000000000e-01,
                                         4.628571428571428571428571428570e+01,
                                         -1.118957142857142857142857142860e+03,
                                         1.067724642857142857142857142860e+04,
                                         -5.356327500000000000000000000000e+04,
                                         1.565208562500000000000000000000e+05,
                                         -2.753159625000000000000000000000e+05,
                                         2.865985392857142857142857142860e+05,
                                         -1.622793053571428571428571428570e+05,
                                         3.843457232142857142857142857140e+04,

                                         0.000000000000000000000000000000e-01,
                                         -1.012500000000000000000000000000e+01,
                                         2.463991071428571428571428571430e+02,
                                         -2.373155357142857142857142857140e+03,
                                         1.205014218750000000000000000000e+04,
                                         -3.574719843750000000000000000000e+04,
                                         6.403125937500000000000000000000e+04,
                                         -6.809087812500000000000000000000e+04,
                                         3.950219933035714285714285714290e+04,
                                         -9.608643080357142857142857142860e+03,

                                         0.000000000000000000000000000000e-01,
                                         1.000000000000000000000000000000e+00,
                                         -2.446071428571428571428571428570e+01,
                                         2.373026785714285714285714285710e+02,
                                         -1.216518750000000000000000000000e+03,
                                         3.652973437500000000000000000000e+03,
                                         -6.643012500000000000000000000000e+03,
                                         7.196596875000000000000000000000e+03,
                                         -4.270508035714285714285714285710e+03,
                                         1.067627008928571428571428571430e+03

};

static double _poly_x_lagrange_coeff9[] = {

    -2.546071428571428571428571428570e+01,
    5.235267857142857142857142857140e+02,
    -4.361464285714285714285714285710e+03,
    1.947796875000000000000000000000e+04,
    -5.147992968750000000000000000000e+04,
    8.303765625000000000000000000000e+04,
    -8.026973437500000000000000000000e+04,
    4.270508035714285714285714285710e+04,
    -9.608643080357142857142857142860e+03,
    0.000000000000000000000000000000e-01,

    8.100000000000000000000000000000e+01,
    -2.666635714285714285714285714290e+03,
    2.760892232142857142857142857140e+04,
    -1.397310750000000000000000000000e+05,
    4.001697421875000000000000000000e+05,
    -6.820159500000000000000000000000e+05,
    6.858910406250000000000000000000e+05,
    -3.758047071428571428571428571430e+05,
    8.647778772321428571428571428570e+04,
    0.000000000000000000000000000000e-01,

    -1.620000000000000000000000000000e+02,
    6.791271428571428571428571428570e+03,
    -8.137592678571428571428571428570e+04,
    4.538207250000000000000000000000e+05,
    -1.391547093750000000000000000000e+06,
    2.493344025000000000000000000000e+06,
    -2.604060900000000000000000000000e+06,
    1.469054764285714285714285714290e+06,
    -3.459111508928571428571428571430e+05,
    0.000000000000000000000000000000e-01,

    2.520000000000000000000000000000e+02,
    -1.132020000000000000000000000000e+04,
    1.469522250000000000000000000000e+05,
    -8.776431000000000000000000000000e+05,
    2.844398531250000000000000000000e+06,
    -5.327696025000000000000000000000e+06,
    5.766134850000000000000000000000e+06,
    -3.348078300000000000000000000000e+06,
    8.071260187500000000000000000000e+05,
    0.000000000000000000000000000000e-01,

    -2.835000000000000000000000000000e+02,
    1.316047500000000000000000000000e+04,
    -1.782131625000000000000000000000e+05,
    1.113994012500000000000000000000e+06,
    -3.769397015625000000000000000000e+06,
    7.336100137500000000000000000000e+06,
    -8.207441943750000000000000000000e+06,
    4.902543225000000000000000000000e+06,
    -1.210689028125000000000000000000e+06,
    0.000000000000000000000000000000e-01,

    2.268000000000000000000000000000e+02,
    -1.073250000000000000000000000000e+04,
    1.491260625000000000000000000000e+05,
    -9.610042500000000000000000000000e+05,
    3.359744578125000000000000000000e+06,
    -6.753729375000000000000000000000e+06,
    7.788932156250000000000000000000e+06,
    -4.782969000000000000000000000000e+06,
    1.210689028125000000000000000000e+06,
    0.000000000000000000000000000000e-01,

    -1.260000000000000000000000000000e+02,
    6.038100000000000000000000000000e+03,
    -8.536083750000000000000000000000e+04,
    5.620043250000000000000000000000e+05,
    -2.014021968750000000000000000000e+06,
    4.158525825000000000000000000000e+06,
    -4.929115275000000000000000000000e+06,
    3.108929850000000000000000000000e+06,
    -8.071260187500000000000000000000e+05,
    0.000000000000000000000000000000e-01,

    4.628571428571428571428571428570e+01,
    -2.237914285714285714285714285710e+03,
    3.203173928571428571428571428570e+04,
    -2.142531000000000000000000000000e+05,
    7.826042812500000000000000000000e+05,
    -1.651895775000000000000000000000e+06,
    2.006189775000000000000000000000e+06,
    -1.298234442857142857142857142860e+06,
    3.459111508928571428571428571430e+05,
    0.000000000000000000000000000000e-01,

    -1.012500000000000000000000000000e+01,
    4.927982142857142857142857142860e+02,
    -7.119466071428571428571428571430e+03,
    4.820056875000000000000000000000e+04,
    -1.787359921875000000000000000000e+05,
    3.841875562500000000000000000000e+05,
    -4.766361468750000000000000000000e+05,
    3.160175946428571428571428571430e+05,
    -8.647778772321428571428571428570e+04,
    0.000000000000000000000000000000e-01,

    1.000000000000000000000000000000e+00,
    -4.892142857142857142857142857140e+01,
    7.119080357142857142857142857140e+02,
    -4.866075000000000000000000000000e+03,
    1.826486718750000000000000000000e+04,
    -3.985807500000000000000000000000e+04,
    5.037617812500000000000000000000e+04,
    -3.416406428571428571428571428570e+04,
    9.608643080357142857142857142860e+03,
    0.000000000000000000000000000000e-01

};

static double _poly_xx_lagrange_coeff9[] = {

    5.235267857142857142857142857140e+02,
    -8.722928571428571428571428571430e+03,
    5.843390625000000000000000000000e+04,
    -2.059197187500000000000000000000e+05,
    4.151882812500000000000000000000e+05,
    -4.816184062500000000000000000000e+05,
    2.989355625000000000000000000000e+05,
    -7.686914464285714285714285714290e+04,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    -2.666635714285714285714285714290e+03,
    5.521784464285714285714285714290e+04,
    -4.191932250000000000000000000000e+05,
    1.600678968750000000000000000000e+06,
    -3.410079750000000000000000000000e+06,
    4.115346243750000000000000000000e+06,
    -2.630632950000000000000000000000e+06,
    6.918223017857142857142857142860e+05,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    6.791271428571428571428571428570e+03,
    -1.627518535714285714285714285710e+05,
    1.361462175000000000000000000000e+06,
    -5.566188375000000000000000000000e+06,
    1.246672012500000000000000000000e+07,
    -1.562436540000000000000000000000e+07,
    1.028338335000000000000000000000e+07,
    -2.767289207142857142857142857140e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    -1.132020000000000000000000000000e+04,
    2.939044500000000000000000000000e+05,
    -2.632929300000000000000000000000e+06,
    1.137759412500000000000000000000e+07,
    -2.663848012500000000000000000000e+07,
    3.459680910000000000000000000000e+07,
    -2.343654810000000000000000000000e+07,
    6.457008150000000000000000000000e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    1.316047500000000000000000000000e+04,
    -3.564263250000000000000000000000e+05,
    3.341982037500000000000000000000e+06,
    -1.507758806250000000000000000000e+07,
    3.668050068750000000000000000000e+07,
    -4.924465166250000000000000000000e+07,
    3.431780257500000000000000000000e+07,
    -9.685512225000000000000000000000e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    -1.073250000000000000000000000000e+04,
    2.982521250000000000000000000000e+05,
    -2.883012750000000000000000000000e+06,
    1.343897831250000000000000000000e+07,
    -3.376864687500000000000000000000e+07,
    4.673359293750000000000000000000e+07,
    -3.348078300000000000000000000000e+07,
    9.685512225000000000000000000000e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    6.038100000000000000000000000000e+03,
    -1.707216750000000000000000000000e+05,
    1.686012975000000000000000000000e+06,
    -8.056087875000000000000000000000e+06,
    2.079262912500000000000000000000e+07,
    -2.957469165000000000000000000000e+07,
    2.176250895000000000000000000000e+07,
    -6.457008150000000000000000000000e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    -2.237914285714285714285714285710e+03,
    6.406347857142857142857142857140e+04,
    -6.427593000000000000000000000000e+05,
    3.130417125000000000000000000000e+06,
    -8.259478875000000000000000000000e+06,
    1.203713865000000000000000000000e+07,
    -9.087641100000000000000000000000e+06,
    2.767289207142857142857142857140e+06,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    4.927982142857142857142857142860e+02,
    -1.423893214285714285714285714290e+04,
    1.446017062500000000000000000000e+05,
    -7.149439687500000000000000000000e+05,
    1.920937781250000000000000000000e+06,
    -2.859816881250000000000000000000e+06,
    2.212123162500000000000000000000e+06,
    -6.918223017857142857142857142860e+05,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01,

    -4.892142857142857142857142857140e+01,
    1.423816071428571428571428571430e+03,
    -1.459822500000000000000000000000e+04,
    7.305946875000000000000000000000e+04,
    -1.992903750000000000000000000000e+05,
    3.022570687500000000000000000000e+05,
    -2.391484500000000000000000000000e+05,
    7.686914464285714285714285714290e+04,
    0.000000000000000000000000000000e-01,
    0.000000000000000000000000000000e-01

};

LagrangeCoeff _lagrange_coeff9(9, _poly_lagrange_coeff9,
                               _poly_x_lagrange_coeff9,
                               _poly_xx_lagrange_coeff9);

} // namespace hiflow
