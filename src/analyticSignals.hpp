#ifndef ANALYTIC_SIGNALS_HPP
#define ANALYTIC_SIGNALS_HPP
#include <stdexcept>
#include <cmath>
namespace
{

/// @brief 
/// @param[in] t    This is the (computed differential) travel time
///                 \f$ T_2 - T_1 \f$ in seconds.
/// @param[in] p1   The first pick time measured in UTC seconds since the epoch.
/// @param[in] p2   The second pick time measured in UTC seconds since the epoch.
/// @param[in] w1   The width of the first boxcar in seconds.
/// @param[in] w2   The width of the second boxcar in seconds.
#pragma omp declare simd uniform(p1, p2, w1, w2)
template<typename T>
[[maybe_unused]]
T analyticBoxcarCorrelation(const T t,
                            const T p1, const T p2,
                            const T w1, const T w2)
{
    T xnorm = 1/(2*(w1*w2));//2*M_PI;
    T dt = p2 - p1;
    T half_w1_p_w2 = 0.5*(w1 + w2);
    T half_w2_m_w1 = 0.5*(w2 - w1);
    T shift1 = t - (dt + half_w1_p_w2);
    T shift2 = t - (dt - half_w1_p_w2);
    T shift3 = t - (dt + half_w2_m_w1);
    T shift4 = t - (dt - half_w2_m_w1);
    T result = shift1*std::copysign(1, shift1) //sgn(shift1) //std::copysign(1, shift1)
             + shift2*std::copysign(1, shift2) //sgn(shift2) //std::copysign(1, shift2)
             - shift3*std::copysign(1, shift3) //sgn(shift3) //std::copysign(1, shift3)
             - shift4*std::copysign(1, shift4); //sgn(shift4);// std::copysign(1, shift4);
    return xnorm*result;
}
template<typename T>
[[maybe_unused]] [[nodiscard]]
T analyticBoxcarCorrelation(const T pick1, const T pick2,
                            const T estimate1, const T estimate2, 
                            const T weight1, const T weight2) 
{
    return ::analyticBoxcarCorrelation(estimate2 - estimate1,
                                       pick1, pick2,
                                       weight1, weight2);
}
/// @param[in] sd1   The standard deviation of the first Gaussian (seconds).
/// @param[in] sd2   The standard devaiation of the second Gaussian (seconds).
/// @result The amplitude normalization factor used in the correlation of two
///         Gaussians:
///         \f [
///             \frac{1}{\sqrt{2 \pi} \sqrt{\sigma_1^2 + \sigma_2^2} }
///         \f ].
template<typename T>
[[maybe_unused]]
T analyticGaussianCorrelationAmplitudeNormalization(const T sd1, const T sd2)
{
    if (sd1 == 0 && sd2 == 0)
    {
        throw std::invalid_argument("Both standard deviations are 0");
    }
    double mag = sd1*sd1 + sd2*sd2;
    double result = 1./std::sqrt(2*M_PI*mag);
    return static_cast<T> (result);
}
/// @param[in] sd1   The standard deviation of the first Gaussian (seconds).
/// @param[in] sd2   The standard devaiation of the second Gaussian (seconds).
/// @result The exponential normalization factor used in the correlation of two
///         Guassians:
///         \f [
///             \frac{1}{2 (\sigma_1^2 + \sigma_2^2) }
///         \f ].
template<typename T>
[[maybe_unused]]
T analyticGaussianCorrelationExponentNormalization(const T sd1, const T sd2)
{
    if (sd1 == 0 && sd2 == 0)
    {
        throw std::invalid_argument("Both standard deviations are 0");
    }
    double result = 1./( 2*(sd1*sd1 + sd2*sd2) );
    return static_cast<T> (result);
}
/// @brief Computes the analytic correlation of two Gaussian pulses.
/// @param[in] t   The (computed differential) time at which to evaluate in
///                seconds, i.e., \f$ T_2 - T_1 \f$.
/// @param[in] p1  The first pick time measured in UTC seconds since the epoch.
/// @param[in] p2  The second pick time measured in UTC seconds since the epoch.
/// @param[in] xnormAmp  The amplitude normalization term, e.g.,:
///                      \f [
///                         \frac{1}{\sqrt{2 \pi} \sqrt{\sigma_1^2 + \sigma_2^2}
///                      \f ].
/// @param[in] xnormExp  The exponential normalization term, e.g.,:
///                      \f [
///                         \frac{1}{2 \sqrt{\sigma_1^2 + \sigma_2^2}
///                      \f ].
/// @result The analytic correlation of two Gaussians evaluated at a time t:
///         \f [
///            g_1(t) \star g_2) (t)
///          = \frac{1}{\sqrt{2 \pi} \sqrt{\sigma_1^2 + \sigma_2^2}
///            \exp \left (
///              -\frac{(t - (t_2 -t_1))^2}{2 \sqrt{\sigma_1^2 + \sigma_2^2}
///            \right )
///         \f ]
#pragma omp declare simd uniform(xnormAmp, xnormExp)
template<typename T>
[[maybe_unused]]
T analyticGaussianCorrelation(const T t,
                              const T p1, const T p2,
                              const T xnormAmp, const T xnormExp)
{
    T dt = p2 - p1; // Differential time
    T res = t - dt;
    return xnormAmp*exp(-(res*res)*xnormExp);
}


}
#endif
