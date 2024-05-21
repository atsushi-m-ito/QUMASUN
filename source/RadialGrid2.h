#pragma once


//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>


/*
動径方向の座標において、グリッド上の半径rの値と、グリッド間隔等に関する情報を保持するクラス

他のクラスに参照させて、メモリは共有で使いまわす
*/
class RadialGrid2{
private:
	double* radius_;
public:
	const int num_grids;
	
	const double xi_min;
	//const double xi_max;
	const double xi_delta;//コンストラクタ初期化の順番の為、必ずxi_min,xi_maxより後に定義すること//
	//const double r_max;
	const double* radius;
	

	/*
	RadialGrid2(const int num_grids_a, const double xi_min_a, const double xi_max_a, const double xi_delta_a) :
		num_grids(num_grids_a),		
		xi_min(xi_min_a),
		//xi_max(xi_max_a),
		xi_delta(xi_delta_a),
		//r_max(exp(xi_max_a)),
		radius_(new double[num_grids_a + 1]),
		radius(radius_)
	{		
	};
	*/

	RadialGrid2(const int num_grids_a, double* radius_assign, const double xi_min_a, const double xi_delta_a) :
		num_grids(num_grids_a),
		xi_min(xi_min_a),
		//xi_max(xi_max_a),
		xi_delta(xi_delta_a),
		//r_max(exp(xi_max_a)),
		radius_(radius_assign),
		radius(radius_)
	{
	};



	~RadialGrid2()
	{
		//delete[] radius_;
	};

	double* GetGridPointer(){
		return radius_;
	}

	/*
	RadialGrid(const int num_grids_a, const double* radius_a) :
		num_grids(num_grids_a), 
		radius(radius_a),
		xi_min(log(radius_a[0])),
		xi_max(log(radius_a[num_grids_a])),
		xi_delta((xi_max - xi_min) / (double)(num_grids_a))
	{};
	*/

	//double r_min;//削除候補//
	//double r_max;//削除候補//

#if 1
	double GetValueByXi(const double xi, const double* values) const {
		return GetValueByXi(xi, values, this->num_grids, this->xi_min, this->xi_delta);
	}

	inline
	static double GetValueByXi(const double xi, const double* values, int num_grids, double xi_min, double xi_delta) {
		const double relative_xi = xi - xi_min;
		const double delta = xi_delta;
		double floor_xi = floor(relative_xi / delta);
		double x = (relative_xi / delta) - floor_xi;
		int i = (int)(floor_xi);

		if (i >= num_grids - 1) {
			return 0.0;
		} else if ((i == 0) || (i == num_grids-2) ){
			const double f0 = values[i];
			const double f1 = values[i + 1];
			
			return f0 * (1.0 - x) + f1 * x;
			
		}

		//xi空間での三次スプライン補完//
		if ((1 < i) && (i < num_grids - 3)) {//二次精度//

			const double fm2 = values[i - 2];
			const double fm1 = values[i - 1];
			const double f0 = values[i];
			const double f1 = values[i + 1];
			const double f2 = values[i + 2];
			const double f3 = values[i + 3];
			
			const double d1f0 = 2.0 / 3.0 * (f1 - fm1) - 1.0 / 12.0 * (f2 - fm2);
			const double d1f1 = 2.0 / 3.0 * (f2 - f0) - 1.0 / 12.0 * (f3 - fm1);
			const double d2f0 = 4.0 / 3.0 * (f1 + fm1) - 1.0 / 12.0 * (f2 + fm2) - 5.0 / 2.0 * f0;
			const double d2f1 = 4.0 / 3.0 * (f2 + f0) - 1.0 / 12.0 * (f3 + fm1) - 5.0 / 2.0 * f1;
			
			auto c3 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return 10.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) - 4.0 * (d1p1 - (d1p0 + d2p0)) + 1.0 / 2.0 * (d2p1 - d2p0);
			};
			auto c4 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return -15.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) + 7.0 * (d1p1 - (d1p0 + d2p0)) - (d2p1 - d2p0);
			};
			auto c5 = [&c3,&c4](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return p1 - p0 - d1p0 - d2p0 / 2.0 - c3(p0, d1p0, d2p0, p1, d1p1, d2p1) - c4(p0, d1p0, d2p0, p1, d1p1, d2p1);
			};
			auto K = [&c3, &c4, &c5](double x, double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return (p0)+x * ((d1p0)+x * ((d2p0 / 2.0) + x * (c3(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * (c4(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * c5(p0, d1p0, d2p0, p1, d1p1, d2p1)))));
			};

			const double x = (relative_xi / delta) - (double)i;
			return K(x, f0, d1f0, d2f0, f1, d1f1, d2f1);

		}
		else /*if ((0 < i) && (i < num_grids - 1))*/ {//二次精度//
			

			const double v0 = values[i];
			const double v1 = values[i + 1];
			const double dvdxi1 = (values[i + 2] - v0) / (2.0);
			const double dvdxi0 = (v1 - values[i - 1]) / (2.0);

			const double a = dvdxi1 + dvdxi0 - 2.0 * (v1 - v0);
			const double b = -dvdxi1 - 2.0 * dvdxi0 + 3.0 * (v1 - v0);
			const double c = dvdxi0;
			const double d = v0;

			const double x = (relative_xi / delta) - (double)i;
			const double x2 = x * x;
			const double x3 = x * x * x;

			const double result = a * x3 + b * x2 + c * x + d;
			/*if(i != (int)floor_xi){
				if(result < 0.0){
					printf("xi=%lg, floor_xi=%lg, i=%d, x=%lg, result=%lg\n", xi, floor_xi, i, x, result);
				}
			}*/
			return result;
			//return v1;
		}
		
		return 0.0;
		
	};

#else
	double GetValueByXi(const double xi, const double* values) const{
		
		const double relative_xi = xi - xi_min;
		const double delta = xi_delta;
		double floor_xi = floor( relative_xi / delta);
		int i = (int)(floor_xi);
		
		if(i==0){ i = 1;}
		if (i == num_grids - 1){ i = num_grids - 2; }

		//xi空間での三次スプライン補完//
		if((0 < i) && (i < num_grids-1)){

			const double v0 = values[i];
			const double v1 = values[i+1];
			const double dvdxi1 = (values[i+2] - v0) / (2.0);
			const double dvdxi0 = (v1 - values[i-1]) / (2.0);
			
			const double a = dvdxi1 + dvdxi0 - 2.0 * (v1-v0);
			const double b = -dvdxi1 - 2.0*dvdxi0 + 3.0 * (v1-v0);
			const double c = dvdxi0;
			const double d = v0;

			const double x = (relative_xi/delta) - (double)i;
			const double x2 = x * x;
			const double x3 = x * x * x;
			
			const double result = a*x3 +b*x2 + c*x + d;
			/*if(i != (int)floor_xi){
				if(result < 0.0){
					printf("xi=%lg, floor_xi=%lg, i=%d, x=%lg, result=%lg\n", xi, floor_xi, i, x, result);
				}
			}*/
			return result;
			//return v1;
		}else{
			return 0.0;
		}	
	};
#endif

	double GetValue(const double r, const double* values) const{
		return GetValueByXi(log(r), values);
	};

	double GetValueBySquare(const double r2, const double* values) const{
		return GetValueByXi(0.5*log(r2), values);
	}

	inline
	static double GetValue(const double r, const double* values, int num_grids, double xi_min, double xi_delta)  {
		return GetValueByXi(log(r), values, num_grids, xi_min, xi_delta);
	};

	inline
	static double GetValueBySquare(const double r2, const double* values, int num_grids, double xi_min, double xi_delta)  {
		return GetValueByXi(0.5 * log(r2), values, num_grids, xi_min, xi_delta);
	}

#if 0	
	double GetValueDebug(const double r, const double* values) const{
		const double xi = log(r);
	

		const double relative_xi = xi - xi_min;
		const double delta = xi_delta;
		double floor_xi = floor( relative_xi / delta);
		double x = (relative_xi / delta) - floor_xi;
		int i = (int)(floor_xi);
		
		if (i == 0) {
			i = 1; x = 0.0;
		}
		if (i == num_grids - 1) {
			i = num_grids - 2;
			x = delta;
		}


		//xi空間での三次スプライン補完//
		if((0 < i) && (i < num_grids-1)){

			const double v0 = values[i];
			const double v1 = values[i+1];
			const double dvdxi1 = (values[i+2] - v0) / (2.0);
			const double dvdxi0 = (v1 - values[i-1]) / (2.0);
			
			const double a = dvdxi1 + dvdxi0 - 2.0 * (v1-v0);
			const double b = -dvdxi1 - 2.0*dvdxi0 + 3.0 * (v1-v0);
			const double c = dvdxi0;
			const double d = v0;

			const double x2 = x * x;
			const double x3 = x * x * x;
			
			const double result = a*x3 +b*x2 + c*x + d;

			//if(i != (int)floor_xi){
			/*	if(result < 0.0){
					printf("r=%lg, xi=%lg, floor_xi=%lg, i=%d, x=%lg, result=%lg\n", r, xi, floor_xi, i, x, result);
					printf("\tvalues=%lg, %lg, %lg, %lg\n", values[i-1], values[i], values[i+1], values[i+2]);
					printf("\ta=%lg, b=%lg, c=%lg, d=%lg\n", a, b, c, d);
				}*/
			//}
			return result;
		}else{
			return 0.0;
		}	
	};
#endif
	/*
	* Pay attention this is not diveded yet.
	*/
	double GetLaplacianByXi(const double xi, const double* values, const int L) const {
		return GetLaplacianByXi(xi, values, L, this->num_grids, this->xi_min, this->xi_delta);

	}

	inline 
	static double GetLaplacianByXi(const double xi, const double* values, const int L, int num_grids, double xi_min, double xi_delta)  {
		const double relative_xi = xi - xi_min;
		const double delta = xi_delta;
		double floor_xi = floor(relative_xi / delta);
		double x = (relative_xi / delta) - floor_xi;
		int i = (int)(floor_xi);
		const double L_L_plus_1 = (double)(L * (L + 1));


		if (i >= num_grids - 2) {
			return 0.0;
		} else if(i <= 1) {
			i = 2; x = 0.0;
		}


		//xi空間での三次スプライン補完//
		if ((1 < i) && (i < num_grids - 3)) {//二次精度//

			const double fm2 = values[i - 2];
			const double fm1 = values[i - 1];
			const double f0 = values[i];
			const double f1 = values[i + 1];
			const double f2 = values[i + 2];
			const double f3 = values[i + 3];

			const double d1f0 = 2.0 / 3.0 * (f1 - fm1) - 1.0 / 12.0 * (f2 - fm2);
			const double d1f1 = 2.0 / 3.0 * (f2 - f0) - 1.0 / 12.0 * (f3 - fm1);
			const double d2f0 = 4.0 / 3.0 * (f1 + fm1) - 1.0 / 12.0 * (f2 + fm2) - 5.0 / 2.0 * f0;
			const double d2f1 = 4.0 / 3.0 * (f2 + f0) - 1.0 / 12.0 * (f3 + fm1) - 5.0 / 2.0 * f1;

			auto c3 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return 10.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) - 4.0 * (d1p1 - (d1p0 + d2p0)) + 1.0 / 2.0 * (d2p1 - d2p0);
			};
			auto c4 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return -15.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) + 7.0 * (d1p1 - (d1p0 + d2p0)) - (d2p1 - d2p0);
			};
			auto c5 = [&c3, &c4](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return p1 - p0 - d1p0 - d2p0 / 2.0 - c3(p0, d1p0, d2p0, p1, d1p1, d2p1) - c4(p0, d1p0, d2p0, p1, d1p1, d2p1);
			};
			auto K = [&c3, &c4, &c5](double x, double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return (p0)+x * ((d1p0)+x * ((d2p0 / 2.0) + x * (c3(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * (c4(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * c5(p0, d1p0, d2p0, p1, d1p1, d2p1)))));
			};
			auto dKdx = [&c3, &c4, &c5](double x, double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return (d1p0)+x * ((d2p0) + x * (3.0*c3(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * (4.0*c4(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * 5.0*c5(p0, d1p0, d2p0, p1, d1p1, d2p1))));
			};
			auto d2Kdx2 = [&c3, &c4, &c5](double x, double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
				return ((d2p0)+x * (6.0 * c3(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * (12.0 * c4(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * 20.0 * c5(p0, d1p0, d2p0, p1, d1p1, d2p1))));
			};

			const double v = K(x, f0, d1f0, d2f0, f1, d1f1, d2f1);
			const double dv_dxi = dKdx(x, f0, d1f0, d2f0, f1, d1f1, d2f1) / delta;
			const double d2v_dxi2 = d2Kdx2(x, f0, d1f0, d2f0, f1, d1f1, d2f1) / (delta * delta);

			double Laplace_Psi = (d2v_dxi2 + dv_dxi - (L_L_plus_1 * v));//
			return -0.5 * Laplace_Psi;


		} else{//二次精度//


			const double v0 = values[i];
			const double v1 = values[i + 1];
			const double dvdxi1 = (values[i + 2] - v0) / (2.0);
			const double dvdxi0 = (v1 - values[i - 1]) / (2.0);

			const double a = dvdxi1 + dvdxi0 - 2.0 * (v1 - v0);
			const double b = -dvdxi1 - 2.0 * dvdxi0 + 3.0 * (v1 - v0);
			const double c = dvdxi0;
			const double d = v0;

			const double x2 = x * x;
			const double x3 = x * x * x;

			const double v = a * x3 + b * x2 + c * x + d;
			
			const double dv_dxi = 3.0*a * x2 + 2.0*b * x + c;
			const double d2v_dxi2 = 6.0 * a * x + 2.0 * b;

			double Laplace_Psi = (d2v_dxi2 + dv_dxi - (L_L_plus_1 * v));//
			return -0.5 * Laplace_Psi;


		} 
			
		return 0.0;
		
	};

	double GetLaplacian(const double r, const double* values, const int L) const {
		return GetLaplacianByXi(log(r), values, L) / (r * r);
	};

	inline
	static double GetLaplacian(const double r, const double* values, const int L, int num_grids, double xi_min, double xi_delta) {
		return GetLaplacianByXi(log(r), values, L, num_grids, xi_min, xi_delta) / (r * r);
	};


    inline
    static double GetDifferential1ByXi(const double xi, const double* values, int num_grids, double xi_min, double xi_delta) {
        const double relative_xi = xi - xi_min;
        const double delta = xi_delta;
        double floor_xi = floor(relative_xi / delta);
        double x = (relative_xi / delta) - floor_xi;
        int i = (int)(floor_xi);
        //const double L_L_plus_1 = (double)(L * (L + 1));


        if (i >= num_grids - 2) {
            return 0.0;
        } else if (i <= 1) {
            i = 2; x = 0.0;
        }


        //xi空間での三次スプライン補完//
        if ((1 < i) && (i < num_grids - 3)) {//二次精度//

            const double fm2 = values[i - 2];
            const double fm1 = values[i - 1];
            const double f0 = values[i];
            const double f1 = values[i + 1];
            const double f2 = values[i + 2];
            const double f3 = values[i + 3];

            const double d1f0 = 2.0 / 3.0 * (f1 - fm1) - 1.0 / 12.0 * (f2 - fm2);
            const double d1f1 = 2.0 / 3.0 * (f2 - f0) - 1.0 / 12.0 * (f3 - fm1);
            const double d2f0 = 4.0 / 3.0 * (f1 + fm1) - 1.0 / 12.0 * (f2 + fm2) - 5.0 / 2.0 * f0;
            const double d2f1 = 4.0 / 3.0 * (f2 + f0) - 1.0 / 12.0 * (f3 + fm1) - 5.0 / 2.0 * f1;

            auto c3 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
                return 10.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) - 4.0 * (d1p1 - (d1p0 + d2p0)) + 1.0 / 2.0 * (d2p1 - d2p0);
                };
            auto c4 = [](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
                return -15.0 * (p1 - (p0 + d1p0 + d2p0 / 2.0)) + 7.0 * (d1p1 - (d1p0 + d2p0)) - (d2p1 - d2p0);
                };
            auto c5 = [&c3, &c4](double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
                return p1 - p0 - d1p0 - d2p0 / 2.0 - c3(p0, d1p0, d2p0, p1, d1p1, d2p1) - c4(p0, d1p0, d2p0, p1, d1p1, d2p1);
                };
            auto dKdx = [&c3, &c4, &c5](double x, double p0, double  d1p0, double  d2p0, double  p1, double  d1p1, double  d2p1) {
                return (d1p0)+x * ((d2p0)+x * (3.0 * c3(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * (4.0 * c4(p0, d1p0, d2p0, p1, d1p1, d2p1) + x * 5.0 * c5(p0, d1p0, d2p0, p1, d1p1, d2p1))));
                };

            const double dv_dxi = dKdx(x, f0, d1f0, d2f0, f1, d1f1, d2f1) / delta;
            
            return dv_dxi;


        } else {//二次精度//


            const double v0 = values[i];
            const double v1 = values[i + 1];
            const double dvdxi1 = (values[i + 2] - v0) / (2.0);
            const double dvdxi0 = (v1 - values[i - 1]) / (2.0);

            const double a = dvdxi1 + dvdxi0 - 2.0 * (v1 - v0);
            const double b = -dvdxi1 - 2.0 * dvdxi0 + 3.0 * (v1 - v0);
            const double c = dvdxi0;
            const double d = v0;

            const double x2 = x * x;
            const double x3 = x * x * x;

            //const double v = a * x3 + b * x2 + c * x + d;

            const double dv_dxi = 3.0 * a * x2 + 2.0 * b * x + c;
            return dv_dxi;


        }

        return 0.0;

    };

    inline
    static double GetDifferential1(const double r, const double* values, int num_grids, double xi_min, double xi_delta) {        
        return GetDifferential1ByXi(log(r), values, num_grids, xi_min, xi_delta) / r ;
        
    };

};

