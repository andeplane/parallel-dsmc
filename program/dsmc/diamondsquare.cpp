#include "diamondSquare.h"
#include "random.h"
#include <math.h>


vector<vector<double > >&DiamondSquare::generate(const uint power2,
        const double H,
        const vector<double> corners,
        const long seed,
        const double sigma,
        const bool addition,
        const bool PBC,
        const uint RNG) {
    this->power2 = power2;
    this->addition = addition;
    this->sigma = sigma;
    this->PBC = PBC;
    this->RNG = RNG;
    systemSize = pow(2.0, power2) + 1;
    rnd = new Random(-(abs(seed)),0,0);
    srand(seed); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()
    R.resize(systemSize);
    for(int i=0; i<systemSize; i++) R[i].resize(systemSize,0);

    if (PBC) { // We need the same value in the corners if we are using periodic boundaries
        R[0][0]                       = corners[0];
        R[0][systemSize-1]            = R[0][0];
        R[systemSize-1][0]            = R[0][0];
        R[systemSize-1][systemSize-1] = R[0][0];
    } else {
        R[0][0]                       = corners[0];
        R[0][systemSize-1]            = corners[1];
        R[systemSize-1][0]            = corners[2];
        R[systemSize-1][systemSize-1] = corners[3];
    }

    runDiamondSquare(R, H, sigma);

    return R;
}

void DiamondSquare::runDiamondSquare(vector<vector<double> >& R, const double H, const double sigma) {

    double RNGstddv = sigma;
    uint stepLength = systemSize-1;
    uint halfStepLength = stepLength/2;

    for (uint depth = 1; depth <= power2; depth++) {

        // Squares
        RNGstddv *= pow(0.5, 0.5*H);
        for (uint x = halfStepLength; x < systemSize - halfStepLength; x += stepLength) {
            for (uint y = halfStepLength; y < systemSize - halfStepLength; y += stepLength) {
                R[x][y] = square(x, y, halfStepLength, RNGstddv, R);
            }
        }

        if (addition) {
            // Add a random number to the points used in the interpolations above
            uint limit;
            if (PBC) {
                limit = systemSize - halfStepLength;
            } else {
                limit = systemSize;
            }
            for (uint x = 0; x < limit; x += stepLength) {
                for (uint y = 0; y < limit; y += stepLength) {
                    R[x][y] += random()*RNGstddv;
                }
            }
            if (PBC) {
                R[0][systemSize-1] = R[0][0];
                R[systemSize-1][0] = R[0][0];
                R[systemSize-1][systemSize-1] = R[0][0];
            }
        }

        // Diamonds
        RNGstddv *= pow(0.5, 0.5*H);
        for (uint x = 0; x < systemSize - halfStepLength; x += stepLength) {
            for (uint y = halfStepLength; y < systemSize - halfStepLength; y += stepLength) {
                R[x][y] = diamond(x, y, halfStepLength, RNGstddv, R);
            }
        }
        for (uint x = halfStepLength; x < systemSize - halfStepLength; x += stepLength) {
            for (uint y = 0; y < systemSize - halfStepLength; y += stepLength) {
                R[x][y] = diamond(x, y, halfStepLength, RNGstddv, R);
            }
        }

        if (PBC) {
            for (uint idx = halfStepLength; idx < systemSize-1; idx+=halfStepLength) {
                R[idx][systemSize-1] = R[idx][0];
                R[systemSize-1][idx] = R[0][idx];
            }
        } else {
            // Bottom edge diamonds
            for (uint y = halfStepLength; y < systemSize - halfStepLength; y += stepLength) {
                uint x = systemSize-1;
                R[x][y] = nonPBCbottomEdgeDiamonds(x, y, halfStepLength, RNGstddv, R);
            }

            // Right edge diamonds
            for (uint x = halfStepLength; x < systemSize - halfStepLength; x+= stepLength) {
                uint y = systemSize-1;
                R[x][y] = nonPBCrightEdgeDiamonds(x, y, halfStepLength, RNGstddv, R);
            }
        }

        if (addition) {
            // Add a random number to the points used in the interpolations above
            uint limit;
            if (PBC) {
                limit = systemSize - halfStepLength;
            } else {
                limit = systemSize;
            }
            for (uint x = 0; x < limit; x += stepLength) {
                for (uint y = 0; y < limit; y += stepLength) {
                    R[x][y] += random()*RNGstddv;
                }
            }
            for (uint x = halfStepLength; x < systemSize-halfStepLength; x += stepLength) {
                for (uint y = halfStepLength; y < systemSize-halfStepLength; y += stepLength) {
                    R[x][y] += random()*RNGstddv;
                }
            }
            if (PBC) {
                R[0][systemSize-1] = R[0][0];
                R[systemSize-1][0] = R[0][0];
                R[systemSize-1][systemSize-1] = R[0][0];
            }
        }

        stepLength /= 2;
        halfStepLength /= 2;
    }
}

double DiamondSquare::square(
        const uint x,
        const uint y,
        const uint halfStepLength,
        const double RNGstddv,
        const vector<vector<double> >&R) {

    return random()*RNGstddv + 0.25*(
        R[x+halfStepLength][y+halfStepLength] +
        R[x+halfStepLength][y-halfStepLength] +
        R[x-halfStepLength][y+halfStepLength] +
        R[x-halfStepLength][y-halfStepLength]);
}

double DiamondSquare::diamond(
        const uint x,
        const uint y,
        const uint halfStepLength,
        const double RNGstddv,
        const vector<vector<double> > &R) {

    double average;

    if (x == 0) { // At top edge of system
        if (PBC) {
            average = 0.25*(
                R[x][y+halfStepLength] +
                R[x][y-halfStepLength] +
                R[x+halfStepLength][y] +
                R[R.size()-1-halfStepLength][y]);
        } else {
            average = (1.0/3.0)*(
                R[x][y+halfStepLength] +
                R[x][y-halfStepLength] +
                R[x+halfStepLength][y]);
        }
    } else if (y == 0) { // At left edge of system
        if (PBC) {
            average = 0.25*(
                R[x][y+halfStepLength] +
                R[x][R[0].size()-1-halfStepLength] +
                R[x+halfStepLength][y] +
                R[x-halfStepLength][y]);
        } else {
            average = (1.0/3.0)*(
                R[x][y+halfStepLength] +
                R[x+halfStepLength][y] +
                R[x-halfStepLength][y]);
        }
    } else {
        average = 0.25*(
            R[x][y+halfStepLength] +
            R[x][y-halfStepLength] +
            R[x+halfStepLength][y] +
            R[x-halfStepLength][y]);
    }

    return average + random()*RNGstddv;
}

double DiamondSquare::nonPBCbottomEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, vector<vector<double> >& R) {

    return random()*RNGstddv + (1.0/3.0)*(
        R[x-halfStepLength][y] +
        R[x][y+halfStepLength] +
        R[x][y-halfStepLength]);
}

double DiamondSquare::nonPBCrightEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, vector<vector<double> >& R) {

    return random()*RNGstddv + (1.0/3.0)*(
        R[x][y-halfStepLength] +
        R[x+halfStepLength][y] +
        R[x-halfStepLength][y]);
}

inline double DiamondSquare::random() {
    // Returns random number with mean 0
    if (RNG == 0) {
        return 0.0;
    } else if (RNG == 1) {
        return rnd->next_double() - 0.5;
        // return (randu<double>() - 0.5); // uniform distribution in [-0.5,0.5]
    } else if (RNG == 2) {
        return rnd->next_gauss();
    } else {
        return NAN;
    }
}
