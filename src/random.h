// random.h - random number generation

#ifndef RANDOM_H
#define RANDOM_H

#include <cpp11.hpp>
#include <random>
#include <functional>
#include <string>

typedef std::mt19937 Engine;

// Discrete class: handles draws from a discrete distribution
class Discrete
{
public:
    Discrete(Engine& engine, cpp11::doubles cdf)
     : e(engine), pmf(cdf.size(), 0), out(cdf.size(), 0)
    {
        // Set up probability mass vector
        for (size_t i = 0; i + 1 < cdf.size(); ++i)
            pmf[i] = cdf[i + 1] - cdf[i];
        pmf[cdf.size() - 1] = 1 - cdf[cdf.size() - 1];
    }

    // Draw new multinomial as conditional binomials
    void operator()(double N)
    {
        out.assign(out.size(), 0);

        double remaining_p = 1.0;
        double remaining_n = N;

        for (size_t k = 0; k < pmf.size() - 1; ++k)
        {
            // Quit if nothing left to assign
            if (remaining_n <= 0 || remaining_p <= 0)
                break;

            // Get next p
            double p = pmf[k] / remaining_p;
            if (p < 0) p = 0;
            if (p > 1) p = 1;

            // Get next bin
            // When variance is above 100, use normal approximation
            double mean = remaining_n * p;
            double var = mean * (1 - p);
            double draw;
            if (var > 100) {
                double x = std::normal_distribution<double>(mean, std::sqrt(var))(e);
                x = std::round(x);
                if (x < 0) x = 0;
                if (x > remaining_n) x = remaining_n;
                draw = x;
            } else {
                draw = (double)std::binomial_distribution<size_t>(remaining_n, p)(e);
            }
            out[k] = draw;
            remaining_n -= draw;
            remaining_p -= pmf[k];
        }

        // Get final bin
        out[out.size() - 1] = remaining_n;
    }

    // Access draw
    double operator[](size_t i) const
    {
        return out[i];
    }

    // Access size of cache
    size_t size() const
    {
        return out.size();
    }

    // Add cached draw to vector; return total number added
    template <class Vec>
    double add_to(Vec& v, size_t start)
    {
        double total = 0;
        for (size_t i = 0; i < out.size() && i + start < v.size(); ++i)
        {
            v[i + start] += out[i];
            total += out[i];
        }
        return total;
    }

private:
    Engine& e;
    std::vector<double> pmf;
    std::vector<double> out;
};

#endif // RANDOM_H
