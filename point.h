#ifndef POINT_H
#define POINT_H


class point
{
public:
    double x;
    double y;
    bool highlighted;

    point(double ax, double ay);
    point(point* pointcopyfrom);
    void setxy(double ax, double ay);
};

#endif // POINT_H
