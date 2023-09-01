#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "body.h"

void xacceleration(body *bodies, int current_body, int n, double G);

void yacceleration(body *bodies, int current_body, int n, double G);

void rungeKutta(body *bodies, int n, double t, double dt, double G);