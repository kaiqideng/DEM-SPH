#include "solverBase.h"

class problem :public solverBase
{
public:
	problem() : solverBase() {}

	void conditionInitialize() override
	{
		// add your code here
	}

	void handleDataAfterCalculateContact() override
	{
		// add your code here
	}

	void outputData() override
	{
		// add your code here
	}
};

class test :public solverBase
{
public:
	test() : solverBase() {}

	void conditionInitialize() override
	{
		setProblemName("test");
		double spacing = 0.015;
		double H = 0.3;
		setDomain(make_double3(-3. * spacing, -3. * spacing, -3. * spacing), make_double3(5.4 * H + 6. * spacing, 2 * H + 6. * spacing, 3 * H + 6. * spacing));
		
		std::vector<double3> fp = getRegularPackedPoints(make_double3(0, 0, 0), make_double3(2 * H, 2 * H, H), spacing);
		addFluid(fp, make_double3(0, 0, 0), 1.3 * spacing, 1000., 10 * 2 * sqrt(9.81 * H), 1.e-6);

		double spacing1 = spacing;
		std::vector<double3> sp0 = getRegularPackedPoints(make_double3(-3 * spacing1, -3 * spacing1, -3 * spacing1), make_double3(5.4 * H + 6 * spacing1, 2 * H + 6 * spacing1, 3 * H + 6 * spacing1), spacing1);
		std::vector<double3> sp1;
		for (const auto& p : sp0)
		{
			if (p.x < 0. || p.y < 0. || p.z < 0. || p.x > 5.4 * H || p.y > 2 * H) sp1.push_back(p);
		}
		addSPHSolidBoundary(sp1, 0.5 * spacing1);
		double dt = 0.25 * 1.3 * spacing / (11 * 2 * sqrt(9.81 * H));
		setSimulationParameterTimeStep(dt);
		setSimulationParameterMaximumTime(5.);
		setSimulationParameterNumFrames(100);
	}

	void outputData() override
	{
		outputFluidVTU();
		outputSolidVTU();
	}
};

class test1 :public solverBase
{
public:
	test1() : solverBase() {}

	void conditionInitialize() override
	{
		setProblemName("test1");
		double spacing = 0.015;
		double H = 0.3;
		setDomain(make_double3(-3. * spacing, -3. * spacing, -3. * spacing), make_double3(5.4 * H + 6. * spacing, 2 * H + 6. * spacing, 3 * H + 6. * spacing));

		std::vector<double3> fp = getRegularPackedPoints(make_double3(0, 0, 0), make_double3(2 * H, 2 * H, H), spacing);
		addFluid(fp, make_double3(0, 0, 0), 1.3 * spacing, 1000., 10 * 2 * sqrt(9.81 * H), 1.e-6);

		double spacing1 = spacing;
		std::vector<double3> sp0 = getRegularPackedPoints(make_double3(- 3 * spacing1, - 3 * spacing1, - 3 * spacing1), make_double3(5.4 * H + 6 * spacing1, 2 * H + 6 * spacing1, 3 * H + 6 * spacing1), spacing1);
		std::vector<double3> sp1;
		for (const auto& p : sp0)
		{
			if (p.x < 0. || p.y < 0. || p.z < 0.|| p.x > 5.4 * H || p.y > 2 * H) sp1.push_back(p);
		}
		addSPHSolidBoundary(sp1, 0.5 * spacing1);
		double dt = 0.25 * 1.3 * spacing / (11 * 2 * sqrt(9.81 * H));
		setSimulationParameterTimeStep(dt);
		setSimulationParameterMaximumTime(5.);
		setSimulationParameterNumFrames(100);

		std::vector<double3> sp2 = getRegularPackedPoints(make_double3(3.4 * H, 0.75 * H, 0), make_double3(0.5 * H, 0.5 * H, 3 * H), spacing1);
		addSPHSolidBoundary(sp2, 0.5 * spacing1);
	}

	void outputData() override
	{
		outputFluidVTU();
		outputSolidVTU();
	}
};

class test2 :public solverBase
{
public:
	test2() : solverBase() {}

	void conditionInitialize() override
	{
		setProblemName("test2");
		double L = 4.;
		int nBond = 10;
		double spacing = L / double(nBond);
		setDomain(make_double3(-L, -0.5 * L, -0.5 * L), make_double3(3 * L , L, L));
		setGravity(make_double3(0, 0, 0));
		setBondedContactModel(0, 0, 200e9, 2.6, 1., 0., 0., 0.);
		std::vector<double3> sp1(nBond + 1);
		for (int i = 0; i < nBond + 1; i++)
		{
			sp1[i] = make_double3(double(i) * spacing, 0, 0);
		}
		std::vector<double> rad(nBond + 1, 0.5 * spacing);
		std::vector<double> den(nBond + 1, 7800.);
		den[0] = 0.;//fixed particle
		addCluster(sp1, rad, den, make_double3(0, 0, 0), 0);
		double dt = 1.e-5;
		setSimulationParameterTimeStep(dt);
		setSimulationParameterMaximumTime(5.);
		setSimulationParameterNumFrames(100);
	}

	void handleDataAfterCalculateContact() override
	{
		addExternalForce(10, make_double3(0, 0, 100.e3));
		for (int i = 0; i <= 10; i++)
		{
			addGlobalDamping(i, 0.1);
		}
	}

	void outputData() override
	{
		const auto& hs = getHostSolidData();
		std::cout << "Z Displacement(m) of each element: " << std::endl;
		std::cout << std::fixed << std::setprecision(6);
		for (int i = 0;i < hs.points.num;i++)
		{
			std::cout << hs.points.position[i].z << " ";
		}
		std::cout << std::endl;
	}
};

class damBreak :public solverBase
{
public:
	damBreak() : solverBase() {}

	double spacing = 0.025;

	void conditionInitialize() override
	{
		setProblemName("damBreak");
		
		setDomain(make_double3(-3. * spacing, -3. * spacing, -3. * spacing), make_double3(8 + 6. * spacing, 0.7 + 6. * spacing, 0.7 + 6. * spacing));

		std::vector<double3> fp = getRegularPackedPoints(make_double3(0, 0, 0), make_double3(3.5, 0.7, 0.4), spacing);
		addFluid(fp, make_double3(0, 0, 0), 1.3 * spacing, 1000., 30., 1.e-6);

		std::vector<double3> sp0 = getRegularPackedPoints(make_double3(-3 * spacing, -3 * spacing, -3 * spacing), make_double3(8 + 6 * spacing, 0.7 + 6 * spacing, 0.8 + 6 * spacing), spacing);
		std::vector<double3> sp1;
		for (const auto& p : sp0)
		{
			if (p.x < 0 || p.y < 0 || p.z < 0 || p.x > 8. || p.y > 0.7) sp1.push_back(p);
		}
		addSolid(sp1, make_double3(0, 0, 0), 0.5 * spacing, 0., 0);

		setHertzianContactModel(0, 1, 3e9, 3e9 / (2 * (1 + 0.3)), 0.9, 0, 0, 0.35, 0, 0);
		setHertzianContactModel(1, 1, 3e9, 3e9 / (2 * (1 + 0.3)), 0.9, 0, 0, 0.45, 0, 0);
		double r_ball = 0.015;
		double k = 0.5 * 3e9 * pi() * r_ball;
		double mij = 4. / 3. * pi() * r_ball * r_ball * r_ball * 800.;
		double dt = 0.2 * sqrt(mij / k);
		double dt_f = 0.25 * 1.3 * spacing / (20 * sqrt(0.4 * 9.8));
		int gap = int(dt_f / dt);
		if (gap < 1) gap = 1;
		setSimulationParameterTimeStep(dt);
		setFluidIntegrateGap(gap);

		setSimulationParameterMaximumTime(2.);
		setSimulationParameterNumFrames(100);
	}

	void handleDataAfterCalculateContact() override
	{
		double r_ball = 0.015;
		double k = 0.5 * 3e9 * pi() * r_ball;
		double mij = 4. / 3. * pi() * r_ball * r_ball * r_ball * 800.;
		double dt = 0.2 * sqrt(mij / k);
		if (getTime() >= 0.5 && getTime() < 0.5 + dt )
		{
			double3 s_cub = make_double3(0.15, 0.15, 0.15);
			std::vector<double3> sp2 = getRegularPackedPoints(make_double3(5.3, 0.275, 0.), s_cub, 2 * r_ball);
			std::vector<double3> sp3 = getRegularPackedPoints(make_double3(5.3, 0.275, 0.15), s_cub, 2 * r_ball);
			std::vector<double3> sp4 = getRegularPackedPoints(make_double3(5.3, 0.275, 0.3), s_cub, 2 * r_ball);
			std::vector<double> rad(sp2.size(), r_ball);
			double mass = 800 * s_cub.x * s_cub.y * s_cub.z;
			symMatrix inertia = make_symMatrix(mass / 12. * (s_cub.y * s_cub.y + s_cub.z * s_cub.z), mass / 12. * (s_cub.x * s_cub.x + s_cub.z * s_cub.z), mass / 12. * (s_cub.x * s_cub.x + s_cub.y * s_cub.y), 0., 0., 0.);
			addClump(sp2, rad, make_double3(5.375, 0.35, 0.075), make_double3(0, 0, 0), mass, inertia, 1);
			addClump(sp3, rad, make_double3(5.375, 0.35, 0.225), make_double3(0, 0, 0), mass, inertia, 1);
			addClump(sp4, rad, make_double3(5.375, 0.35, 0.375), make_double3(0, 0, 0), mass, inertia, 1);
		}
	}

	void outputData() override
	{
		outputFluidVTU();
		outputSolidVTU();
	}
};

int main()
{
	test1 problem;
	problem.solve();
}