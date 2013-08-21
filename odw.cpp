/*
 *
 * To use this program, copy it to Test_Skies/odw.cpp and
 * run the following command:
 *
 * make "CPPFLAGS=-O3" odw && ./odw > results.csv
 *
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

// The model used has two parameters. Optimum values were
// found using cross validation on the training data.
const double ALPHA  = 0.22;
const double BETA   = 125;

double sqr(double val)
{
    return (val * val);
}

class Point
{
public:
    Point() : x(0), y(0)
    {
    }

    void print()
    {
        printf(",%.2f,%.2f", x, y);
    }

    double x, y;
};

class Galaxy
{
public:
    void read(FILE* file)
    {
        int id;
        int result = fscanf(file, "Galaxy%d,%lf,%lf,%lf,%lf\n",
                            &id, &x, &y, &e1, &e2);  
        assert(result == 5);
    }

    double distance(Point& p)
    {
        return sqrt(sqr(p.x - x) + sqr(p.y - y)); 
    }


    double  x, y, e1, e2;
};

class Sky
{
public:
    Sky(int id)
    {
        char    buffer[1024];
        sprintf(buffer, "Test_Sky%d.csv", id+1);
        FILE*   file = fopen(buffer, "r");
        assert(file != 0);
        fscanf(file, "%s\n", buffer);

        for (galaxyCount = 0; feof(file) == 0; galaxyCount++) {
            assert(galaxyCount < sizeof(galaxies) / sizeof(galaxies[0]));
            galaxies[galaxyCount].read(file);
        }

        haloCount = (id < 40) ? 1 : ((id < 80) ? 2 : 3); 
        fclose(file);
    }

    void findHalos(Point* halos)
    {
        for (unsigned int k = 0; k < haloCount; k++) {
            if (k > 0) {
                wipeHalo(halos[k-1]);
            }

            findNextHalo(halos[k]);
        }
    }

private:
    double getSignalStrength(Point& p)
    {
        double  sum = 0;
        double  weightSum = 0;

        for (unsigned int i = 0; i < galaxyCount; i++) {
            Galaxy& g = galaxies[i];
            double  r = g.distance(p);
            if (r >= 4200) {
                continue;
            }

            double  weight = 4200 - r;
            double  phi = atan2(g.y - p.y, g.x - p.x);
            sum -= weight * (g.e1 * cos(2 * phi) + g.e2 * sin(2 * phi)); 
            weightSum += weight;
        }

        // Return the weighted average of tangential
        // ellipticity observed from the given point.
        return (sum / weightSum);
    }

    void wipeHalo(Point& h)
    {
        // Estimate the mass of the halo to be a multiple
        // of the signal strength seen at its center.
        double  mass = getSignalStrength(h) * BETA;
        for (unsigned int i = 0; i < galaxyCount; i++) {
            Galaxy& g = galaxies[i];
            double  phi = atan2(g.y - h.y, g.x - h.x);
            double  r = g.distance(h);
            // Model the gravitational pull according to the halo's estimated
            // mass and the distance from the center of the halo. 
            double  force = mass * exp(-pow(r, ALPHA));
            g.e1 += force * cos(2 * phi);
            g.e2 += force * sin(2 * phi);
        }
    }

    void findNextHalo(Point& h)
    {
        Point   ref;
        int     sizes[] = {105, 21, 3, 1};
        int     curr = 4200;

        // Divide the sky into square tiles and check the signal strength
        // at the tile's center. Pick the tile with the strongest
        // signal, divide it further and iterate.
        for (unsigned int iter = 0; iter < 4; iter++) { 
            double  max = -1;
            int     count = curr / sizes[iter];  
            curr = sizes[iter];

            for (int row = 0; row < count; row++) {
                for (int col = 0; col < count; col++) {
                    Point   p;
                    p.x = ref.x + col * curr + curr / 2;
                    p.y = ref.y + row * curr + curr / 2;
                    double  sig = getSignalStrength(p);
                    if (sig > max) { 
                        max = sig; 
                        h = p;
                    }
                }
            }

            ref.x = h.x - curr / 2;
            ref.y = h.y - curr / 2;
        }
    }

    Galaxy          galaxies[800];
    unsigned int    galaxyCount;
    unsigned int    haloCount;
};

int main()
{
    for (int j = 0; j < 120; j++) {
        Point   halos[3];
        Sky     sky(j);

        sky.findHalos(halos);
        printf("Sky%d", j+1);
        for (int k = 0; k < 3; k++) {
            halos[k].print();
        }

        printf("\n");
    }

    return 0;
}

