//
// Created by xetql on 1/23/19.
//

#ifndef CA_ROAD_METRIC_HPP
#define CA_ROAD_METRIC_HPP

#include <utility>

template<typename Realtype, typename ContainerA, typename ContainerB>
std::pair<Realtype, Realtype> linear_regression(const ContainerA& x, const ContainerB& y) {
    int i; Realtype xsomme, ysomme, xysomme, xxsomme;

    Realtype ai, bi;

    xsomme = 0.0; ysomme = 0.0;
    xysomme = 0.0; xxsomme = 0.0;
    const int n = x.size();
    for (i=0;i<n;i++) {
        xsomme = xsomme + x[i]; ysomme = ysomme + y[i];
        xysomme = xysomme + x[i]*y[i];
        xxsomme = xxsomme + x[i]*x[i];
    }
    ai = (n*xysomme - xsomme*ysomme)/(n*xxsomme - xsomme*xsomme);
    bi = (ysomme - ai*xsomme)/n;

    return std::make_pair(ai, bi);
}

template<typename Realtype, typename ContainerA, typename ContainerB>
Realtype get_slope(const ContainerA& x, const ContainerB& y){

    int n = x.size();

    Realtype avgX = std::accumulate(x.begin(), x.end(), 0.0) / n;
    Realtype avgY = std::accumulate(y.begin(), y.end(), 0.0) / n;

    Realtype numerator = 0.0;
    Realtype denominator = 0.0;

    for(int i=0; i<n; ++i){
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    return numerator / denominator;
}
#endif //CA_ROAD_METRIC_HPP
