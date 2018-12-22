//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_DRIVING_DIRECTION_HPP
#define CA_ROAD_DRIVING_DIRECTION_HPP

enum DrivingDirection {
    GoingLeft, GoingRight, GoingUp, GoingDown, Rotary, NoDirection
}; // with int: 0, 1, 2, 3, 4

DrivingDirection char_to_driving(char c) {
    switch (c) {
        case 'L':
            return GoingLeft;
        case 'R':
            return GoingRight;
        case 'U':
            return GoingUp;
        case 'D':
            return GoingDown;
        case 'O':
            return Rotary;
        default :
            return NoDirection;
    }
}

char driving_to_char(DrivingDirection c) {
    switch (c) {
        case GoingLeft:
            return 'L';
        case GoingRight:
            return 'R';
        case GoingUp:
            return 'U';
        case GoingDown:
            return 'D';
        case Rotary:
            return 'O';
        default :
            return '-';
    }
}

#endif //CA_ROAD_DRIVING_DIRECTION_HPP
