#pragma once
#include "integrate.h"
#include "fileEditingTool.h"
#include "generatePoints.h"

struct simulationParameter
{
    double maximumTime;
    double timeStep;
    int numSteps;
    int iStep;
    int numFrames;
    int iFrame;;
    int frameInterval;
};

struct GPUParameter
{
    int deviceIndex;
    int maxThreadsPerBlock;
};

class data
{
public:
    data()
    {
        simPara.maximumTime = 1.;
        simPara.timeStep = 1.;
        simPara.numSteps = 1;
        simPara.iStep = 0;
        simPara.numFrames = 1;
        simPara.iFrame = 0;
        simPara.frameInterval = 1;
        gpuPara.deviceIndex = 0;
        gpuPara.maxThreadsPerBlock = 256;
        domainOrigin = make_double3(0, 0, 0);
        domainSize = make_double3(1, 1, 1);
        gravity = make_double3(0, 0, -9.81);
    }

    ~data()
    {
        dev.release();
    }

protected:

    void setGPUParameterDeviceIndex(int index)
    {
        gpuPara.deviceIndex = index;
    }

    void setSimulationParameterMaximumTime(double t)
    {
        simPara.maximumTime = t;
    }

    void setSimulationParameterTimeStep(double dt)
    {
        simPara.timeStep = dt;
    }

    void setSimulationParameterNumFrames(int n)
    {
        simPara.numFrames = n;
    }

    void setDomain(double3 origin, double3 size)
    {
        domainOrigin = origin;
        domainSize = size;
        if (simPara.iStep > 1) setSpatialGrids();
    }

    void setGravity(double3 g)
    {
        gravity = g;
        if (simPara.iStep > 1) dev.gravity = gravity;
    }

    void setHertzianContactModel(size_t mat_i, size_t mat_j, double E, double G, double res, double k_r_k_s, double k_t_k_s, double mu_s, double mu_r, double mu_t)
    {
		int c_ij = hos.contactModels.getCombinedIndex(mat_i, mat_j);
		if (c_ij < 0)
		{
			std::cout << "Error: material index exceeds the number of materials." << std::endl;
			return;
		}
		hos.contactModels.hertzian.E[c_ij] = E;
		hos.contactModels.hertzian.G[c_ij] = G;
		hos.contactModels.hertzian.res[c_ij] = res;
		hos.contactModels.hertzian.k_r_k_s[c_ij] = k_r_k_s;
		hos.contactModels.hertzian.k_t_k_s[c_ij] = k_t_k_s;
		hos.contactModels.hertzian.mu_s[c_ij] = mu_s;
		hos.contactModels.hertzian.mu_r[c_ij] = mu_r;
		hos.contactModels.hertzian.mu_t[c_ij] = mu_t;
    }

	void setLinearContactModel(size_t mat_i, size_t mat_j, double k_n, double k_s, double k_r, double k_t, double d_n, double d_s, double d_r, double d_t, double mu_s, double mu_r, double mu_t)
	{
		int c_ij = hos.contactModels.getCombinedIndex(mat_i, mat_j);
        if (c_ij < 0)
        {
            std::cout << "Error: material index exceeds the number of materials." << std::endl;
            return;
        }
		hos.contactModels.linear.k_n[c_ij] = k_n;
		hos.contactModels.linear.k_s[c_ij] = k_s;
		hos.contactModels.linear.k_r[c_ij] = k_r;
		hos.contactModels.linear.k_t[c_ij] = k_t;
		hos.contactModels.linear.d_n[c_ij] = d_n;
		hos.contactModels.linear.d_s[c_ij] = d_s;
		hos.contactModels.linear.d_r[c_ij] = d_r;
		hos.contactModels.linear.d_t[c_ij] = d_t;
		hos.contactModels.linear.mu_s[c_ij] = mu_s;
		hos.contactModels.linear.mu_r[c_ij] = mu_r;
		hos.contactModels.linear.mu_t[c_ij] = mu_t;
	}

	void setBondedContactModel(size_t mat_i, size_t mat_j, double E, double k_n_k_s, double gamma, double sigma_s, double C, double mu)
	{
		int c_ij = hos.contactModels.getCombinedIndex(mat_i, mat_j);
        if (c_ij < 0)
        {
            std::cout << "Error: material index exceeds the number of materials." << std::endl;
            return;
        }
		hos.contactModels.bonded.E[c_ij] = E;
		hos.contactModels.bonded.k_n_k_s[c_ij] = k_n_k_s;
		hos.contactModels.bonded.gamma[c_ij] = gamma;
		hos.contactModels.bonded.sigma_s[c_ij] = sigma_s;
		hos.contactModels.bonded.C[c_ij] = C;
		hos.contactModels.bonded.mu[c_ij] = mu;
	}

    void addFluid(std::vector<double3> p, double3 velocity, double smoothLength, double density, double soundSpeed, double kinematicViscosity);

    void addSolid(std::vector<double3> p, double3 velocity, double radius, double density, int materialID);

	void addCluster(std::vector<double3> p, std::vector<double3> velocity, std::vector<double> radius, std::vector<double> density, int materialID);

    void addClump(std::vector<double3> p, std::vector<double> radius, double3 centroidPosition, double3 velocity, double mass, symMatrix inertiaTensor, int materialID);

    void addExternalForce(int index, double3 force);

    void addGlobalDamping(int index, double C_d);

	const double getTime() const
	{
		return simPara.iStep * simPara.timeStep;
	}

	const int getStep() const
	{
		return simPara.iStep;
	}

	const int getFrame() const
	{
		return simPara.iFrame;
	}

private:
    HostData hos;
    DeviceData dev;
    simulationParameter simPara;
    GPUParameter gpuPara;
    double3 domainOrigin;
    double3 domainSize;
    double3 gravity;

    void addFluidData(const HostFluid f);

    void addSolidData(const HostSolid s);

    void setSpatialGrids();

    void buildDeviceData();

	friend class solverBase;
};