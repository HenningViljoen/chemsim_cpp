#ifndef POINT_H
#define POINT_H


class point
{
public:
    double x;
    double y;
    bool highlighted;

    point();
    point(double ax, double ay);
    point(point const &pointcopyfrom);

    point &operator=(const point &p);

    void copyfrom(const point &pointcopyfrom);
    void setxy(double ax, double ay);
};

#endif // POINT_H
