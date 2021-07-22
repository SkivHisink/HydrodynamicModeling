#pragma once
#include "CombustionModel.hpp"


class FireLine
{
public:
	//Начальные значения переменных
	double P = 1e10;//10GPa
	double cv = 1300;// Дж/(кг*К)
	double heat = 4e6; // Теплота реакции Дж/кг
	double lambda = 10;// теплопроводность Вт/(м*К)
	//double ro1 = 1.8e3; //средняя плотность до нагружения
	double ro2 = 2e3; //средняя плотность после нагружения
	double E = 17900;// Температура активации в К
	double time_ = 0;
	double k_ = 1e3;//
//	double k;
	double alpha = 1;
	const double P_ = 1e5;
	double tau = 1e-10;//
	double dx = 1e-6;
	double timeStep_mult = 0.00005;
	double addToTimeStep = 0.005;
	int lineSize = 1000;
	int front = 0;
	//
	std::vector<double> Temp;
	std::vector<double> N;
	std::vector<double> N_mean, time_N;
	std::vector<double> coordinates;
	std::vector<double> dN_val;
	std::vector<double> Temp_mean;
	std::vector<double> combustion_time;
	double max(std::vector<double> val)
	{
		double max_val = val[0];
		for (int i = 0; i < val.size(); ++i)
		{
			if (max_val < val[i]) { max_val = val[i]; }
		}
		return max_val;
	}
	double max_tau(const double max) const { return 0.1 / (k_ * (1 + alpha * P / P_)) * exp(E / max); }
	double mean(const std::vector<double>& arr)
	{
		double summ = 0.0;
#pragma omp parallel for collapse(1) reduction(+:summ)
		for (int k = 0; k < arr.size(); ++k)
		{
			summ += arr[k];
		}
		return summ / arr.size();
	}
	int substance_exist(const std::vector<double>& N)
	{
		auto mean_ = mean(N);
		if (mean_ < 0.01)
			return 0;
		return 1;
	}
	void init(const std::string& density_source, int rand_size)
	{
		//amount of unreacted substance in point
		N.assign(lineSize, 1.0);
		//Start temperature in Kelvin
		Temp.assign(lineSize, 300);
		for (int i = 0; i < 30; ++i)
			Temp[i] = 4000;

		N_mean.push_back(1);
		time_N.push_back(0);
		dN_val.push_back(0);
		Temp_mean.push_back(300);
		for (int i = 0; i < lineSize; ++i)
		{
			coordinates.push_back(i);
		}

	}
	void simulate(const std::string& temp_name, const std::string& N_name)
	{
		double time_period = time_;
		bool scheme_state = true;
		std::vector<double> Temp_new = Temp;
		const double temp_tau = 1.0 / 2 * pow(dx, 2) * cv * ro2 / lambda;
		int sub = 1;
		double mean_val = mean(N);
		int step = 0;
		double old_time = time_;
		int prev_front = front;
		while (substance_exist(N))
		{
			double tmp = 0.0;
			// вычисление шага по времени
			const double max_tau_val = max_tau(max(Temp));
			if (max_tau_val < temp_tau)
				tau = max_tau_val;
			else
				tau = temp_tau;

#pragma omp parallel for collapse(1) reduction(+:tmp)
			for (int k = 1; k < Temp.size() - 1; ++k)
			{
				double dN = -N[k] * k_ * (1 + alpha * P / P_) * exp(-E / Temp[k]) * tau;
				if (isnan(dN) || isinf(dN))
				{
					scheme_state = false; //i = j = k = INT32_MAX - 2;
					continue;
				}
				Temp_new[k] = Temp[k] - heat / cv * dN +
					tau * lambda / (pow(dx, 2) * cv * ro2) *
					(Temp[k - 1] + Temp[k + 1] - 2 * Temp[k]);
				if (isnan(Temp_new[k]) || isinf(Temp_new[k]))
				{
					scheme_state = false; //i = j = k = INT32_MAX - 2;
					continue;
				}
				N[k] += dN;
				tmp += dN;
			}
			if (!scheme_state)
			{
				break;
			}
			//std::cout << "N_mean: " << mean_val << " |Time: " << time_ << " |tau: " << tau << std::endl;
			Temp = Temp_new;

			Temp[0] = Temp[1];
			N[0] = N[1];
			Temp[Temp.size() - 1] = Temp[Temp.size() - 2];

			// измерение скорости волны горения
			if (step % 1000 == 0)
			{
				for (int i = 1; i < N.size(); ++i)
				{
					if (N[i] > 0.5)
					{
						front = i;
						break;
					}
				}
				std::cout << "Speed: " << (double)(front - prev_front) * dx / (time_ - old_time) << std::endl;
				//
				/*combustion_time.push_back(time_);
				timeStep_mult += addToTimeStep;
				std::string save_path = "E:\\University\\6 semester\\Institute\\fireline\\07202021\\";
				plt::scatter(coordinates, Temp);
				plt::show();
				plt::save(save_path + "Temp" + std::to_string((double)N.size() / sub / lineSize) + ".png");
				plt::scatter(coordinates, N);
				plt::save(save_path + "N" + std::to_string((double)N.size() / sub / lineSize) + ".png");
				plt::show();*/
				//
				prev_front = front;
				old_time = time_;
			}
			if (N[N.size() / 8 * sub] <= 0.5)
			{
				if (sub != 9)
				{
					++sub;
					combustion_time.push_back(time_);
					timeStep_mult += addToTimeStep;
					std::string save_path = "E:\\University\\6 semester\\Institute\\fireline\\07202021\\";
#ifndef _DEBUG
					plt::scatter(coordinates, Temp);
					plt::show();
					plt::save(save_path + "Temp" + std::to_string((double)N.size() / sub / lineSize) + ".png");
					plt::scatter(coordinates, N);
					plt::save(save_path + "N" + std::to_string((double)N.size() / sub / lineSize) + ".png");
					plt::show();
#endif
				}
			}

			time_ += tau;
			if (time_ >= 1000e-6 * timeStep_mult)
			{
				//
			}
			if (time_period < time_)
			{
				N_mean.push_back(mean_val);
				time_N.push_back(time_);
				time_period = time_;
				dN_val.push_back(tmp / N.size());
				Temp_mean.push_back(mean(Temp));
			}
			mean_val = mean(N);
			++step;
		}
#ifndef _DEBUG
		//std::string save_path = "E:\\University\\6 semester\\Institute\\fireline\\07202021\\";
		//plt::scatter(coordinates, Temp);
		//plt::save(save_path + "Temp_final" + ".png");
		//plt::show();
		//plt::scatter(coordinates, N);
		//plt::show();
#endif
		for (int i = 0; i < combustion_time.size() - 1; ++i)
		{
			std::cout << combustion_time[i] << " ";
		}
		std::cout << time_;
	}

};