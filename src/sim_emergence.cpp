#include <cpp11.hpp>
#include "random.h"

[[cpp11::register]]
cpp11::data_frame simulate_emergence(cpp11::doubles R, double k,
    cpp11::doubles generation_cdf,
    cpp11::list outcome_cdfs,
    unsigned int seed)
{
    // Simulation length is length of R(t) vector
    size_t time = R.size();

    // Create RNG engine and discrete distributions
    Engine e(seed);
    Discrete infection_dist(e, generation_cdf);
    std::vector<Discrete> outcome_dists;
    for (size_t i = 0; i < outcome_cdfs.size(); ++i)
        outcome_dists.emplace_back(e, outcome_cdfs[i]);

    // Number of secondary infections with R = r and n source infections
    auto secondary_dist = [&](double r, double n) -> double {
        if (n <= 0.0 || r <= 0.0) return 0.0;

        std::gamma_distribution<double> g(k * n, 1.0 / k);
        double lambda = r * g(e);
        if (!std::isfinite(lambda))
            cpp11::stop("Non-finite lambda in secondary_dist().");

        // Above lambda = 1e6, use normal approximation
        if (lambda < 1e6)
            return std::poisson_distribution<size_t>(lambda)(e);
        double x = std::normal_distribution<double>(lambda, std::sqrt(lambda))(e);
        return std::max(0.0, std::round(x));
    };

    // Create infections time series
    auto infections = cpp11::writable::doubles(time);
    std::fill(infections.begin(), infections.end(), 0);
    infections[0] = 1;
    double active_infections = 1;

    // Create outcomes time series
    std::vector<cpp11::writable::doubles> outcomes;
    for (size_t i = 0; i < outcome_cdfs.size(); ++i)
    {
        outcomes.emplace_back(time);
        std::fill(outcomes.back().begin(), outcomes.back().end(), 0);
    }

    // Loop through each time step
    for (size_t t = 0; t < time; ++t)
    {
        // cpp11::check_user_interrupt();

        // Add secondary infections
        double nsec = secondary_dist(R[t], infections[t]);
        active_infections -= infections[t];
        infection_dist(nsec);
        active_infections += infection_dist.add_to(infections, t);

        // Account for secondary infections added "today"
        size_t extra = 0;
        while (infection_dist[0] > 0) {
            nsec = secondary_dist(R[t], infection_dist[0]);
            active_infections -= infection_dist[0];
            infection_dist(nsec);
            active_infections += infection_dist.add_to(infections, t);

            ++extra;
            if (extra > 20)
                cpp11::stop("Supercritical reproduction.");
        }

        // Add outcomes
        for (size_t i = 0; i < outcomes.size(); ++i)
        {
            outcome_dists[i](infections[t]);
            outcome_dists[i].add_to(outcomes[i], t);
        }

        // Quit early if no more active infections
        if (active_infections <= 0)
            break;
    }

    // Return infections and outcomes as a data.frame
    cpp11::writable::list result;
    result.push_back(infections);
    for (size_t i = 0; i < outcomes.size(); ++i)
        result.push_back(outcomes[i]);

    // Apply names: "infections" followed by names from outcome_cdfs
    cpp11::writable::strings names;
    names.push_back("infections");
    SEXP onames = outcome_cdfs.attr("names");
    if (onames == R_NilValue) {
        for (size_t i = 0; i < outcomes.size(); ++i)
            names.push_back("");
    } else {
        cpp11::strings on(onames);
        for (size_t i = 0; i < outcomes.size(); ++i)
            names.push_back(on[i]);
    }
    result.attr("names") = names;

    // Promote to data.frame (class + compact row.names)
    result.attr("class") = "data.frame";
    cpp11::writable::integers row_names({NA_INTEGER, -static_cast<int>(time)});
    result.attr("row.names") = row_names;

    return cpp11::data_frame(result);
}
