//
// Created by xetql on 12/21/18.
//

#ifndef CA_ROAD_TCA_UTILS_HPP
#define CA_ROAD_TCA_UTILS_HPP

#include <map>
#include <unordered_map>
#include <algorithm>

template<class key, class val>
val get_or_default(const std::unordered_map<key, val> map, const key k1, const val def){
    if(map.find(k1) == map.cend())
        return def;
    return map.at(k1);
}
inline long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}
inline long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}
inline std::pair<int, int> cell_to_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msy);
}

#endif //CA_ROAD_TCA_UTILS_HPP
