//
// Created by Алексей on 19.11.2024.
//
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <fstream>

class Point {
public:
  double x_;
  double y_;

  Point(double x, double y) : x_(x), y_(y) {};
  double squareDistance(Point p) {
    return std::pow(p.x_ - x_, 2) + std::pow(p.y_ - y_, 2);
  }
  double distance(Point p) {return sqrt(squareDistance(p)); }
  double X() { return x_; }
  double Y() { return y_; }
};

class Circle {
private:
  Point point_;
  double radius_;
public:
  Circle(double x, double y, double r) : point_(x, y), radius_(r) {};
  bool pointInCircle(Point p) {
    return p.squareDistance(point_) <= pow(radius_, 2);
  }

  double X() {return point_.X(); }
  double Y() {return point_.Y(); }
  double R() {return radius_; }
};

void bigRectangleBoundaries(Circle first, Circle second, Circle third,
  double &leftBound, double &rightBound, double &lowerBound, double &upperBound) {

  leftBound = std::min(first.X() - first.R(),
    std::min(second.X() - second.R(), third.X() - third.R()));
  rightBound = std::max(first.X() + first.R(),
  std::min(second.X() + second.R(), third.X() + third.R()));
  lowerBound = std::min(first.Y() - first.R(),
    std::min(second.Y() - second.R(), third.Y() - third.R()));
  upperBound = std::max(first.Y() + first.R(),
  std::min(second.Y() + second.R(), third.Y() + third.R()));
}


 // Сгенерируем точку внутри области и от нее построим прямоугольник.
void smallRectangleBoundaries(Circle first, Circle second, Circle third, int seed,
  double &leftBound, double &rightBound, double &lowerBound, double &upperBound) {
  bigRectangleBoundaries(first, second, third, leftBound,
    rightBound, lowerBound, upperBound);
  leftBound = std::max(first.X() - first.R(),
    std::min(second.X() - second.R(), third.X() - third.R()));
  rightBound = std::min(first.X() + first.R(),
  std::min(second.X() + second.R(), third.X() + third.R()));
  lowerBound = std::max(first.Y() - first.R(),
    std::min(second.Y() - second.R(), third.Y() - third.R()));
  upperBound = std::min(first.Y() + first.R(),
  std::min(second.Y() + second.R(), third.Y() + third.R()));
}

double estimateArea(Circle first, Circle second, Circle third,
  int numberOfTests, int seed, double leftBound,
  double rightBound, double lowerBound, double upperBound) {

  std::mt19937 mt(seed);
  int count = 0;
  for (int i = 0; i < numberOfTests; ++i) {

    double x = leftBound + (mt() * 1.0 / std::mt19937::max()) * (rightBound - leftBound);
    double y = lowerBound + (mt() * 1.0  / std::mt19937::max()) * (upperBound - lowerBound);

    Point p = Point(x, y);
    if (first.pointInCircle(p) && second.pointInCircle(p) &&
      third.pointInCircle(p)) {
      ++count;
      }
  }

  return count / (numberOfTests * 1.0) *
    ((rightBound - leftBound) * (upperBound - lowerBound));
}

double estimatedArea_bigRectangle(Circle first, Circle second, Circle third,
  int numberOfTests, int seed) {
  double leftBound, rightBound, lowerBound, upperBound;
  bigRectangleBoundaries(first, second, third, leftBound,
    rightBound, lowerBound, upperBound);

  return estimateArea(first, second, third, numberOfTests, seed, leftBound,
    rightBound, lowerBound, upperBound);
}

double estimatedArea_smallRectangle(Circle first, Circle second, Circle third,
  int numberOfTests, int seed) {
  double leftBound, rightBound, lowerBound, upperBound;
  smallRectangleBoundaries(first, second, third, seed, leftBound,
    rightBound, lowerBound, upperBound);

  return estimateArea(first, second, third, numberOfTests, seed, leftBound,
    rightBound, lowerBound, upperBound);
}


int mainA() {
  // double x, y, r;
  // std::cin >> x >> y >> r;
  // auto first = Circle(x, y, r);
  // std::cin >> x >> y >> r;
  // auto second = Circle(x, y, r);
  // std::cin >> x >> y >> r;
  // auto third = Circle(x, y, r);

  double x = 1, y = 1, r = 1;
  auto first = Circle(x, y, r);
  x = 1.5, y = 2, r = sqrt(5)/2;
  auto second = Circle(x, y, r);
  x = 2, y = 1.5, r = sqrt(5)/2;
  auto third = Circle(x, y, r);
  double answer = 0.9445171858994636001513014244729117261206349855590771147710;

  std::ofstream outSmall, outBig, smallRelation, bigRelation;
  outSmall.open("smallRectangle1000000.txt");
  outBig.open("bigRectangle1000000.txt");
  smallRelation.open("smallRelation1000000.txt");
  bigRelation.open("bigRelation1000000.txt");
  outSmall << "x;y\n";
  outBig << "x;y\n";
  smallRelation << "x;y\n";
  bigRelation << "x;y\n";

  for (int i = 100; i <= 100000; i += 500) {
    double small = estimatedArea_smallRectangle(first, second, third, i, 1000000);
    double big = estimatedArea_bigRectangle(first, second, third, i, 1000000);
    outSmall << i << ';' << small << '\n';
    outBig << i << ';' << big << '\n';
    smallRelation << i << ';' << 100 * std::abs(small - answer) / answer << '\n';
    bigRelation << i << ';' << 100 * std::abs(big - answer) / answer << '\n';
  }
  outSmall.close();
  outBig.close();
  smallRelation.close();
  bigRelation.close();
}