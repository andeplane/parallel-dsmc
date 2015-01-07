#ifndef DIAMONDSQUARE_H
#define DIAMONDSQUARE_H

#include <iostream>
#include <vector>
using namespace std;
typedef unsigned int uint;
class Random;

/*! \brief This class creates a randomly generated heightmap using either successive random additions algorithm or
           successive random displacement algorithm.

    The program starts with generating a 2d grid of size (2^power2+1)x(2^power2+1). We initialize the system by setting
    the z-value in the corners to a given value. We then generate the rest of the z-values by iteration.

    The algorithm consists of two separate steps, called "diamond" and "square".

    In the diamond-step we find the z-value of a point (x,y) in the grid by taking the average of four points:
    (x+d, y), (x-d,y), (x,y+d), (x,y-d), where d is a steplength, and adding a random number (usually a normally
    distributed random number, but uniform distribution is also available). The four points make up a diamond centered
    around (x,y), thus we call it the diamond-step. See the ASCII-illustration below.

        o o D o o
        o o o o o
        D o X o D
        o o o o o
        o o D o o
      Figure: The diamond-step. We interpolate the points marked "D" and add a random number to find the value in "X"

    In the square-step we find the z-value of a point (x,y) average of four diagonal neighbors of the point: (x+d, y+d),
    (x+d,y-d), (x-d,y+d), (x-d,y-d). These four neighbors make up a square witht he point we want to find in the center,
    thus we call it the square-step. See the ASCII-illustration below.

        S o o o S
        o o o o o
        o o X o o
        o o o o o
        S o o o S
      Figure: The square-step. We interpolate the points marked "S" and add a random number to find the value in "X"

    When generating the heightmap we first do the square step once, interpolating the four initial points. We then do the
    diamond-step four times, filling in the values in the center of the four edges of the grid. We then continue doing
    this until the whole grid is filled in. The standard deviation of the random displacement in each step is multiplied
    with sigma_(i-1)*0.5^(0.5*H), where sigma_(i-1) is the standard deviation in the last step. This is to produce
    Brownian motion with a Hurst-exponent of H.

    By default the program also adds a random displacement to the points used in the interpolation. This is called the
    successive random additions algorithm. For more info on what this does see the book "The Science of Fractal Images",
    chapter 2 "Algorithms for random fractals", section 2.3.3 (available here:
    http://bringerp.free.fr/Files/Captain%20Blood/Saupe87d.pdf). To turn of this behaviour, and use the successive
    random displacement algorithm instead, set the boolean input variable "addition" to false.

    The methods used to produce the heightmaps are more thoroughly described in the book "Fractals" by Jens Feder,
    sections 9.8 and 13.4 (1988 version), in the book "The Science of Fractal Images", chapters 1 and 2 (1988 version),
    and in the book "Fundamental Algorithms for Computer Graphics" in the chapter "Random fractal forgeries" by Voss R. F.
*/

class DiamondSquare {
public:
    /*!
        \brief generate The method used to generate a heightmap. Returns a reference to the member armadillo matrix.

        This method generates the heightmap as described in the main description of this class.

        The size of the grid is calculated as 2^\a power2 + 1.

        If periodic boundaries is enabled the values in the bottom and top row are the same, as are the values in the
        right and left column. The biggest impact this has is that when we do the diamond step for a point close to or
        at the edge, we will use a point from the opposite end of the grid in the interpolation, instead of just using
        three points. See the ASCII-art below for a visual explanation

            1 o o o # o o o #
            o o o o o o o o o
            D o 2 o o o 4 o o
            o o o o o o o o o
            3 o o o # o o o #
            o o o o o o o o o
            o o # o o o # o o
            o o o o o o o o o
            # o o o # o o o #
        diamond-square    with PBC: the point "D" is calculated by interpolating the points 1,2,3, and 4.
        diamond-square without PBC: the point "D" is calculated by interpolating the points 1,2, and 3.
        (when using PBC all four corners of the grid would be point 1, and points 3 and D would also appear on the right
        edge)

        If addition is enabled we add a random displacement to all the points used in the interpolation for the new
        points in both the diamond- and the square-step. This is called the successive random addition algorithm. If
        addition is disabled we use the successive random displacement algorithm.

        If periodic boundary conditions are enabled, the first element in \a corners is used for all corners.

        \param power2 an integer argument that sets the size of the grid.
        \param H the Hurst-exponent
        \param corners an armadillo vector with the initial values of the corners.
        \param seed the seed used for the random number generator.
        \param sigma the initial standard deviation of the random displacement.
        \param addition a bool controls random displacement of the points used in the interpolation of a new point.
               If true we use the successive random additions algorithm, if false we use the successive random
               displacement algorithm.
        \param PBC a bool that controls periodic boundary conditions.
        \param RNG an unsigned in the selects which random number generator to use (0 just returns 0.0, 1 uses a uniform
               distribution, 2 uses a normal distribution).
    */
    vector<vector<double> >& generate(const uint power2,
            const double H,
            const vector<double> corners,
            const long seed,
            const double sigma,
            const bool addition,
            const bool PBC,
            const uint RNG);

private:
    /*!
        The function that does the diamond- and square-steps.
    */
    void runDiamondSquare(vector<vector<double> > &R, const double H, const double sigma);
    double square(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, const vector<vector<double> > &R);
    double diamond(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, const vector<vector<double> > &R);

    /*!
        Since we can toggle the periodic boundaries we have made separate methods to do the diamond-step without
        periodic boundaries. In the main loop we do the inner points, and the top and left edge. If we are using PBC the
        bottom and right edge have the same values as the top and left edge, so we don't have to do any calculations on
        the bottom and right edge. If we don't have PBC we have to calculate these edges outside the main loop. To
        optimize this a bit we have made separate methods for this.

        \param x the first coordinate of the point we want to calculate.
        \param y the second coordinate of the point we want to calculate.
        \param halfStepLength half the steplength used in this iteration/at this depth.
        \param RNGstddv the standard deviation of the random displacement.
        \param R the armadillo matrix that has the current heightmap.
        \return returns the value of the new point.
     */
    double nonPBCbottomEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, vector<vector<double> > &R);
    /*!
        \brief rightEdgeDiamonds See \a nonPBCbottomEdgeDiamonds.
        \sa nonPBCbottomEdgeDiamonds()
    */
    double nonPBCrightEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, vector<vector<double> > &R);

    /*!
        \brief random the function that generates the random number for the random displacement.
        For now we're just using the built-in methods of Armadillo.
        \return returns a random number, whose distribution depends on the boolean member variable \a PBC.
     */
    double random();

    vector<vector<double> > R;
    uint power2, systemSize, zerolength;
    int RNG;
    bool PBC;
    double sigma;
    bool addition;
    Random *rnd;
};

#endif // DIAMONDSQUARE_H
