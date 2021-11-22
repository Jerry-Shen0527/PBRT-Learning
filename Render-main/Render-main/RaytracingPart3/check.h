#pragma once
#include <iostream>

static class check {
public:
	static void CHECK(bool event);

	static void DCHECK(bool event);

	static void CHECK_LT(double v1, double v2);

	static void CHECK_GT(double v1, double v2);

	static void CHECK_GE(double v1, double v2);

	template <typename T, typename U, typename V>
	static T Clamp(T val, U low, V high);

	static double Lerp(double t, double v1, double v2);
};
