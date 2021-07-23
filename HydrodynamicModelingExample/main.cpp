#include "CombustionModel.hpp"
#include "FireLine.hpp"

enum typeOfModeling
{
	CombustionParallelepiped,
	CombustionCylinder,
	Fireline
};

int main(int argc, char* argv[])
{
	int typeOfM = CombustionParallelepiped;
	std::string destanation = "rhotnt_p_clip.vtk";
	if (typeOfM == CombustionParallelepiped || typeOfM == CombustionCylinder)
	{
		CombustionModel* model = new CombustionModel();
		array3d<unsigned short> density;
		model->initDensity(destanation, density);
		model->init(density, typeOfM);
		model->setSaveDataStep(true);
		model->setShowTimeEachEpoch(true);
		model->simulate("temp.txt", "N.txt");
		model->showReactionDuration();
		model->drawPlotN();
		model->drawPlotdN();
		model->drawPlotTemp();
	}
	if (typeOfM == Fireline) {
		FireLine* line = new FireLine();
		line->init("", 0);
		line->simulate("", "");
	}
	return 0;
}

////std::string source = "rhotnt_p_cylinder.vtk";
	//std::string destanation = "rhotnt_p_clip.vtk";
	//array3d<unsigned short> temp;
	//array3d<unsigned short> temp2(50, 50, 50);
	///*vtk2array3d(source, temp);
	//for (int i = 0; i < temp2.x();++i)
	//for (int j = 0; j < temp2.y();++j)
	//for (int k = 0; k < temp2.z();++k)
	//{
	//	temp2(i, j, k) = temp(temp.x() / 2 + i, temp.y() / 2 + j, temp2.z() / 2 + k);
	//}
	//array3d2vtk(temp2, destanation);*/
	//CombustionModel* model = new CombustionModel();
	//model->init(destanation);
	//model->setSaveDataStep(true);
	//model->setShowTimeEachEpoch(true);
	//model->simulate("temp.txt", "N.txt");
	//model->showReactionDuration();
	//model->drawPlotN();
	//model->drawPlotdN();
	//model->drawPlotTemp();
