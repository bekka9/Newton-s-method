#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> 
using namespace std;
long double f(long double x) {
    return(x * x - sin(10 * x));
}
long double df(long double x) {
    return(2 * x - 10 * cos(10 * x));
}

long double f1(long double x, long double y, long double lambda) {
    return (lambda*cos(y + 0.5) - 2 - x);
}
long double f2(long double x, long double y, long double lambda) {
    return (lambda*sin(x) - 1 - 2 * y);
}
long double f1dx(long double x, long double y, long double lambda) {
    return (-1);
}
long double f2dx(long double x, long double y, long double lambda) {
    return (lambda*cos(x));
}
long double f1dy(long double x, long double y, long double lambda) {
    return (-lambda*sin(y + 0.5));
}
long double f2dy(long double x, long double y, long double lambda) {
    return (-2);
}

void newton_system() {
    long double x = 0, y = 0, lambda = 0, n = 1, h, g, eps = 0.00001;
    x = lambda * cos(y + 0.5) - 2;
    y = (lambda * sin(x) - 1) / 2;
    //F1 = cos(y + 0.5) - 2 - x;
    //F2 = sin(x) - 1 - 2*y;
    long double delta, x_new, y_new;    
    for (int i = 1; i <= n; i++) {
        lambda = i/n;
        h = (f2(x, y, lambda) * f1dx(x, y, lambda) - f1(x, y, lambda) * f2dx(x, y, lambda)) / (f1dy(x, y, lambda) * f2dx(x, y, lambda) - f1dx(x, y, lambda) * f2dy(x, y, lambda));
        g = -(f1(x, y, lambda) + f1dy(x, y, lambda) * h) / f1dx(x, y, lambda);
        x = x + g;
        y = y + h;
    }
    delta = 1; 
    while (delta >= eps) {
        h = (f2(x, y, lambda) * f1dx(x, y, lambda) - f1(x, y, lambda) * f2dx(x, y, lambda)) / (f1dy(x, y, lambda) * f2dx(x, y, lambda) - f1dx(x, y, lambda) * f2dy(x, y, lambda));
        g = -(f1(x, y, lambda) + f1dy(x, y, lambda) * h) / f1dx(x, y, lambda);
        x_new = x + g;
        y_new = y + h;
        delta = sqrt((x_new - x) * (x_new - x) + (y_new - y) * (y_new - y));
        //delta = max(abs(x_new - x), abs(y_new - y));
        x = x_new;
        y = y_new;
    }
    cout << x << ' ' << y;
}
void newton_system_graphic() {
    long double x = -1, y = -1, lambda = 1, h, g, eps = 0.00001, x_new, y_new;
    long double delta = 1;
    while (delta >= eps) {
        h = (f2(x, y, lambda) * f1dx(x, y, lambda) - f1(x, y, lambda) * f2dx(x, y, lambda)) / (f1dy(x, y, lambda) * f2dx(x, y, lambda) - f1dx(x, y, lambda) * f2dy(x, y, lambda));
        g = -(f1(x, y, lambda) + f1dy(x, y, lambda) * h) / f1dx(x, y, lambda);
        x_new = x + g;
        y_new = y + h;
        delta = sqrt((x_new - x) * (x_new - x) + (y_new - y) * (y_new - y));
        
        x = x_new;
        y = y_new;
    }
    cout << x << ' ' << y;
}
void newton_euation() {
    long double a = -5, b = 4, first_otr = 0, x = a, g, eps = 0.0001, x_new = 0;
    //последовательный перебор для определения начального интервала
    int n = 10;
    while (first_otr == 0) {
        for (int i = 0; i < n; i++) {
            x_new = a + i * (b - a) / n;
            if (f(x) * f(x_new) < 0) {
                a = x; b = x_new;
                first_otr = 1;
            }
            else x = x_new;
        }
        n *= 2;
    }

    cout << "x_0 = " << x << '\n';

    //метод ньютона в связке с методом половинного деления для уточнения корня на выбранном интервале локализации
    long double delta = 1;
    x_new = x - f(x) / df(x);
    if ((x_new >= a) && (x_new <= b)) {
        if (f(x_new) < 0) a = f(x_new);
        else b = f(x_new);
    }
    else {
        cout << "x_new = " << x_new << "  a = " << a << "  b = " << b << '\n';
        x_new = (a + b) / 2;
        if (f(x_new) < 0) a = f(x_new);
        else b = f(x_new);
       
       // a = -5, b = 4, n = 10;
        /*-0.68
        - 0.584862 - 0.68 - 0.59
        - 0.592454*/
    }
    while (delta >= eps) {
        x_new = x - f(x) / df(x);
        delta = abs(x_new - x);
        x = x_new;

    }
    cout << "ans = " << x;
}
int main()
{
 /*   newton_system();
    cout << '\n';
    newton_system_graphic();*/
    newton_euation();
}
