#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <random>
#include <string>
#include <omp.h>

#include "array3d.hpp"
#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

class CombustionModel
{
protected:
	enum typeOfModeling
	{
		Parallelepiped,
		Cylinder
	};
	//Начальные значения переменных
	double P = 1e10;//10GPa
	double cv = 1300;// Дж/(кг*К)
	double heat = 4e6; // Теплота реакции Дж/кг
	double lambda = 10;// теплопроводность Вт/(м*К)
	double ro1 = 1.6e3; //средняя плотность до нагружения
	double ro2 = 2.4e3; //средняя плотность после нагружения
	double E = 17900;// Температура активации в К
	double time_ = 0;
	double k_B = 1e11;//
	double tau = 1e-10;//
	double dx = 1e-8;
	//Данные
	array3d<double> density;
	array3d<double> N;
	array3d<double> Temp;
	//Для графиков
	std::vector<double> N_mean, time_N;
	std::vector<double> dN_val;
	std::vector<double> Temp_mean;
	//Вспомогательные переменные
	bool show_logs = true;
	bool model_init = false;
	bool saveMeanVal = true;
	bool saveStepData = true;
	bool showTimeEachEpoch = true;
	double timeStep_mult = 0.00005;
	double addToTimeStep = 0.005;//
	int typeOfFigure = Parallelepiped;
public:
	CombustionModel() = default;
	virtual void init(const std::string& density_source = "", int rand_size = 10);
	void saveInitDensity(const std::string& destination) const { array3d2vtk(density, destination); }
	void saveInitTemperature(const std::string& destination) const { array3d2vtk(Temp, destination); }
	virtual void simulate(const std::string& temp_name = "", const std::string& N_name = "");
	void showReactionDuration() const { std::cout << "Duration of reaction: " << time_ << std::endl; }
	void drawPlotN()
	{
#ifndef _DEBUG
		plt::plot(time_N, N_mean);
		plt::show();
#endif
	}
	void drawPlotdN()
	{
#ifndef _DEBUG
		plt::scatter(time_N, dN_val);
		plt::show();
#endif
	}
	void drawPlotTemp()
	{
#ifndef _DEBUG
		plt::scatter(time_N, Temp_mean);
		plt::show();
#endif
	}
	//set values
	void setP(const double P_) { P = P_; }
	void setCv(const double cv_) { cv = cv_; }
	void setHeat(const double heat_) { heat = heat_; }
	void setLambda(const double lambda_) { lambda = lambda_; }
	void setRo1(const double ro1_) { ro1 = ro1_; }
	void setRo2(const double ro2_) { ro2 = ro2_; }
	void setE(const double E_) { E = E_; }
	void setK_B(const double k_B_) { k_B = k_B_; }
	void setDx(const double dx_) { dx = dx_; }
	void setShowLogs(const bool show_) { show_logs = show_; }
	void setSaveMeanVal(const bool saveMeanVal_) { saveMeanVal = saveMeanVal_; }
	void setSaveDataStep(const bool saveDataStep_) { saveStepData = saveDataStep_; }
	void setAddToTimeStep(const double addToTimeStep_) { addToTimeStep = addToTimeStep_; }
	void setShowTimeEachEpoch(const bool showTimeEachEpoch_) { showTimeEachEpoch = showTimeEachEpoch_; }
	//get values
	double getP() const { return P; }
	double getCv() const { return cv; }
	double getHeat() const { return heat; }
	double getLambda() const { return lambda; }
	double getRo1() const { return ro1; }
	double getRo2() const { return ro2; }
	double getE() const { return E; }
	double getTime() const { return time_; }
	double getK_B() const { return k_B; }
	double getTau() const { return tau; }
	double getDx() const { return dx; }
protected:
	double max_tau(const double max) const { return 0.1 / k_B * exp(E / max); }
	void save_data(const array3d<double>& arr, std::string& name) const;
};

namespace
{
	const double random_in_range(double left_boundary_, double right_boundary_, std::mt19937 rand_engine)
	{
		std::uniform_real_distribution<double> dist(left_boundary_, right_boundary_);
		return dist(rand_engine);
	}
	template<class T>
	bool init_density(array3d<T>& density, int size_, double ro1, const std::string& source = "")
	{
		//if type==0 then create random density from 0.9 to 1.1
		if (source.empty()) {
			std::random_device rand_device;
			//array3d<T> temp(size_, size_, size_);
			array3d<T> temp(3, 3, 500);
			for (int i = 0; i < temp.size(); ++i)
			{
				temp(i) = random_in_range(0.9 * ro1, ro1, std::mt19937(rand_device()));
			}
			density = temp;
			return true;
		}
		return vtk2array3d(source, density);
	}

	bool init_density(array3d<int>& density, int size_, double ro1, const std::string& source = "")
	{
		//if source.empty() then create random density from 0.9 to 1.1
		if (source.empty()) {
			std::random_device rand_device;
			//array3d<T> temp(size_, size_, size_);
			array3d<int> temp(size_, size_, size_);
#pragma omp parallel for collapse(1)
			for (int i = 0; i < temp.size(); ++i)
			{
				temp(i) = 2000/* random_in_range(0.9 * ro1, ro1, std::mt19937(rand_device()))*/;
			}
			density = temp;
			return true;
		}
		return vtk2array3d(source, density);
	}

	double mean(const array3d<double>& arr)
	{
		double summ = 0.0;
#pragma omp parallel for collapse(3) reduction(+:summ)
		for (int k = 1; k < arr.z() - 1; ++k)
		{
			for (int j = 1; j < arr.y() - 1; ++j)
			{
				for (int i = 1; i < arr.x() - 1; ++i)
				{
					summ += arr(i, j, k);
				}
			}
		}
		return summ / ((arr.x() - 2) * (arr.y() - 2) * (arr.z() - 2));
	}

	int substance_exist(const array3d<double>& N)
	{
		auto mean_ = mean(N);
		if (mean_ < 0.01)
			return 0;
		return 1;
	}

}

inline void CombustionModel::init(const std::string& density_source, int typeOfFigure_)
{
	typeOfFigure = typeOfFigure_;
	if (show_logs) { std::cout << "Starting of model initialization" << std::endl; }
	//density data
	if (show_logs) { std::cout << "!!WARNING!! Now density use <int> parametr" << std::endl; }
	array3d<unsigned short> density_tmp;
	//initialization of density
	if (show_logs) { std::cout << "Initialization of density..." << std::endl; }
	if (!init_density(density_tmp, 10, ro1, density_source))
	{
		if (show_logs) { std::cout << "Problem with density initialization." << std::endl; }
		return;
	}
	if (show_logs) { std::cout << "Density successful initialized." << std::endl; }
	density = array3d<double>(density_tmp.x(), density_tmp.y(), density_tmp.z());
	double sum = 0;
#pragma omp parallel for collapse(1) reduction(+:sum)
	for (int i = 0; i < density_tmp.size(); ++i)
	{
		sum += (double)density_tmp(i);
	}
	double normal = sum / density_tmp.size();
#pragma omp parallel for collapse(1)
	for (int i = 0; i < density_tmp.size(); ++i)
	{
		density(i) = ((double)density_tmp(i) + 1.0) / (normal)*ro1;
	}
	//amount of unreacted substance in point
	N = array3d<double>(density.x(), density.y(), density.z(), 1.0);
	//Start temperature in Kelvin
	Temp = array3d<double>(density.x(), density.y(), density.z());
#pragma omp parallel for collapse(1)
	for (int i = 0; i < density.size(); ++i)
	{
		double temp_val = 300 + P / cv * (1 / density(i) - 1 / ro2);
		if (temp_val > 5000.0) { Temp(i) = 5000.0; }
		else if (temp_val < 270.0) { Temp(i) = 300.0; }
		else { Temp(i) = temp_val; }
	}
	N_mean.push_back(mean(N));
	time_N.push_back(0);
	dN_val.push_back(0);
	array3d<double> new_Temp = array3d<double>(Temp.x(), Temp.y(), Temp.z() / 1.5);
	for (int i = 0; i < new_Temp.x(); ++i)
	{
		for (int j = 0; j < new_Temp.y(); ++j) {
			for (int k = 0; k < new_Temp.z(); k += 2) {
				if (3 * k / 2 + 3 < Temp.z())
				{
					auto max_tmp = std::max(Temp(i, j, 3 * k / 2 + 1), Temp(i, j, 3 * k / 2 + 2));
					new_Temp(i, j, k) = Temp(i, j, 3 * k / 2);
					new_Temp(i, j, k + 1) = max_tmp;
				}
				else if (3 * k / 2 + 2 < Temp.z())
				{
					new_Temp(i, j, k) = Temp(i, j, 3 * k / 2);
					new_Temp(i, j, k + 1) = Temp(i, j, 3 * k / 2 + 1);
				}
				else if (3 * k / 2 + 1 < Temp.z())
				{
					new_Temp(i, j, k) = Temp(i, j, 3 * k / 2);
				}
			}
		}
	}
	array3d2vtk(new_Temp, "temp_init" + std::to_string(P) + ".vtk");
	array3d2vtk(Temp, "temp2_init" + std::to_string(P) + ".vtk");
	array3d2vtk(density, "density_init" + std::to_string(P) + ".vtk");
	Temp = new_Temp;
	Temp_mean.push_back(mean(Temp));
	array3d2vtk(Temp, "temp_init.vtk");
	tau = max_tau(Temp.max()); //ограничение на шаг по времени с химической реакцией
	model_init = true;
	if (show_logs) { std::cout << "Model successful initialized." << std::endl; }
}

inline void CombustionModel::simulate(const std::string& temp_name, const std::string& N_name)
{
	double time_period = time_;
	bool scheme_state = true;
	array3d<double> Temp_new(Temp);
	const double temp_tau = 1.0 / 6 * pow(dx, 2) * cv * ro2 / lambda;
	int sub = 8;
	double mean_val = mean(N);
	bool figureCondition = true;
	double radius = Temp.x() / 2;
	while (substance_exist(N))
	{
		double tmp = 0.0;
#pragma omp parallel for collapse(3) reduction(+:tmp)
		for (int k = 1; k < Temp.z() - 1; ++k)
		{
			for (int j = 1; j < Temp.y() - 1; ++j)
			{
				for (int i = 1; i < Temp.x() - 1; ++i)
				{
					if (typeOfFigure == Cylinder)
					{
						figureCondition = sqrt(pow(radius - i, 2) + pow(radius - j, 2)) < radius;
					}
					if (figureCondition) {
						double dN = -N(i, j, k) * k_B * exp(-E / Temp(i, j, k)) * tau;
						if (isnan(dN) || isinf(dN))
						{
							if (show_logs)std::cout << "The scheme is destroyed. " << std::endl << "dN is:" << dN <<
								std::endl << "Temp is:" << Temp(i, j, k) << std::endl;
							scheme_state = false; //i = j = k = INT32_MAX - 2;
							continue;
						}
						Temp_new(i, j, k) = Temp(i, j, k) - heat / cv * dN +
							tau * lambda / (pow(dx, 2) * cv * ro2) *
							(Temp(i - 1, j, k) + Temp(i + 1, j, k) + Temp(i, j - 1, k) + Temp(i, j + 1, k) +
								Temp(i, j, k - 1) + Temp(i, j, k + 1) - 6 * Temp(i, j, k));
						if (isnan(Temp_new(i, j, k)) || isinf(Temp_new(i, j, k)))
						{
							if (show_logs)std::cout << "The scheme is destroyed. " << std::endl << "dN is:" << dN <<
								std::endl << "Temp is:" << Temp(i, j, k) << std::endl;
							scheme_state = false; //i = j = k = INT32_MAX - 2;
							continue;
						}
						N(i, j, k) += dN;
						tmp += dN;
					}
				}
			}
		}
		if (!scheme_state)
		{
			break;
		}
		if (showTimeEachEpoch)
			std::cout << "N_mean: " << mean_val << " |Time: " << time_ << std::endl;
		Temp = Temp_new;
		const double max_tau_val = max_tau(Temp.max());
		if (max_tau_val < temp_tau)
			tau = max_tau_val;
		else
			tau = temp_tau;

		if (saveStepData)
		{
			if (time_ >= 1e-6 * timeStep_mult)
			{
				std::string time_str = temp_name + std::to_string(timeStep_mult) + "_mult_1e-6";
				save_data(Temp, time_str);
				time_str = N_name + std::to_string(timeStep_mult) + "_mult_1e-6";
				save_data(N, time_str);
				timeStep_mult += addToTimeStep;
			}
		}
		time_ += tau;
		for (int j = 0; j < Temp.y(); ++j)
		{
			for (int i = 0; i < Temp.x(); ++i)
			{
				Temp(i, j, 0) = Temp(i, j, 1);
				Temp(i, j, Temp.z() - 1) = Temp(i, j, Temp.z() - 2);
			}
		}
		for (int k = 0; k < Temp.z(); ++k)
		{
			for (int i = 0; i < Temp.x(); ++i)
			{
				Temp(0, i, k) = Temp(1, i, k);
				Temp(Temp.x() - 1, i, k) = Temp(Temp.x() - 2, i, k);
				Temp(i, 0, k) = Temp(i, 1, k);
				Temp(i, Temp.y() - 1, k) = Temp(i, Temp.y() - 2, k);
			}
		}
		Temp(0, 0, 0) = Temp(0, 0, 1);
		Temp(0, 0, Temp.z() - 1) = Temp(0, 0, 1);
		Temp(0, Temp.y() - 1, 0) = Temp(0, 0, 1);
		Temp(0, Temp.y() - 1, Temp.z() - 1) = Temp(0, 0, 1);
		Temp(Temp.x() - 1, 0, 0) = Temp(0, 0, 1);
		Temp(Temp.x() - 1, 0, Temp.z() - 1) = Temp(0, 0, 1);
		Temp(Temp.x() - 1, Temp.y() - 1, 0) = Temp(0, 0, 1);
		Temp(Temp.x() - 1, Temp.y() - 1, Temp.z() - 1) = Temp(0, 0, 1);
		array3d2vtk(Temp, "problem.vtk");
		if (saveMeanVal)
		{
			if (time_period < time_)
			{

				N_mean.push_back(mean_val);
				time_N.push_back(time_);
				time_period = time_;
				dN_val.push_back(tmp / N.size());
				Temp_mean.push_back(mean(Temp));
			}
		}
		mean_val = mean(N);
	}
}

inline void CombustionModel::save_data(const array3d<double>& arr, std::string& name) const
{
	array3d<double> A(arr.x() - 2, arr.y() - 2, arr.z() - 2);
	for (size_t k = 1; k < A.z() + 1; ++k)
		for (size_t j = 1; j < A.y() + 1; ++j)
			for (size_t i = 1; i < A.x() + 1; ++i)
			{
				A(i - 1, j - 1, k - 1) = arr(i, j, k);
			}
	std::string file_name = name + ".vtk";
	if (array3d2vtk(A, file_name.c_str()))
	{
		if (show_logs) std::cout << "Data saved correctly" << std::endl;
	}
	else
	{
		if (show_logs) std::cout << "Something went wrong in data saving" << std::endl;
	}
}
