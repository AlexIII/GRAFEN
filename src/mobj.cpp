#include "mobj.h"

void toRad(limits &l) {
	l.lower = toRad(l.lower);
	l.upper = toRad(l.upper);
}
double toRad(const double a) {
	return M_PI*(a / 180.);
}
void toDeg(limits &l) {
	l.lower = toDeg(l.lower);
	l.upper = toDeg(l.upper);
}
double toDeg(const double a) {
	return 180.*(a / M_PI);
}

std::ostream& operator<<(std::ostream& os, const limits& l) {
	os << "{" << l.lower << ", " << l.upper << ", " << l.n << "}";
	return os;
}
std::ostream& operator<<(std::ostream& os, const Point& p) {
	os << "{" << p.x << ", " << p.y << ", " << p.z << "}";
	return os;
}
std::ostream& operator<<(std::ostream& os, const Triangle& q) {
	os << q.p1 << ", " << q.p2 << ", " << q.p3;
	return os;
}
std::ostream& operator<<(std::ostream& os, const Quadrangle& q) {
	os << q.p1 << ", " << q.p2 << ", " << q.p3 << ", " << q.p4;
	return os;
}
std::ostream& operator<<(std::ostream& os, const PointValue& p) {
	os << p.val << " " << p.x << " " << p.y << " " << p.z;
	return os;
}
std::istream& operator>>(std::istream& is, const PointValue& p) {
	is >> p.val >> p.x >> p.y >> p.z;
	return is;
}