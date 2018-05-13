#include "Search.h"

PointerTwoDim Search::Two_Dim_Search(const double _Epsilon, const int _Steps, const int threads, vagrish::GrishaginFunction *func)
{

	PointerTwoDim p;
	p.steps = 0;
	p.time = 0;

	int procs = threads;

	int localSteps = 0; // шаги для внутреннего счета

	double eps = _Epsilon;

	//текущие границы  единичного отрезка
	double curr_left;
	double curr_right;

	// время работы
	double start, end;

	double *x_new = new double[procs];

	double *coordinates_y = new double[2];

	//константы для оценки конст Липшеца
	double m_small = 1;
	double M_big;

	std::vector<GroupTwoDim> StartVec(procs);
	std::vector<GroupTwoDim> MainVec(procs * 2);
	std::list<GroupTwoDim> Local_group;

	size_t size_group; //// размер list<Group_2D> Local_group

	std::list<GroupTwoDim>::iterator Pointer;

	double  right_x_new, left_x_new;

	omp_set_dynamic(0);//среда выполнения не будет динамически настраивать количество потоков


	std::vector<TEvolvent*> evolve(procs);//вектор разверток для потоков (у каждого своя)

	for (int i = 0; i < procs; i++)//инициализация разверток каждому потоку
	{
		evolve[i] = new TEvolvent(2, 20);
		evolve[i]->SetBounds(Left, Right);
	}

	left_x_new = x_new[0] = 0;
	right_x_new = x_new[procs - 1] = 1;

	double startCooord[2];


	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)//распределение отрезка по потокам
		{

			StartVec[i].x_new_left = left_x_new + (i * (right_x_new - left_x_new) / procs);
			StartVec[i].x_new_right = StartVec[i].x_new_left + (right_x_new - left_x_new) / procs;

			evolve[i]->GetImage(StartVec[i].x_new_left, startCooord);

			StartVec[i].z_left = Func(startCooord, func);

			evolve[i]->GetImage(StartVec[i].x_new_right, startCooord);

			StartVec[i].z_right = Func(startCooord, func);
		}

		// подсчет M
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			StartVec[i].M = M(StartVec[i].z_right, StartVec[i].z_left, StartVec[i].x_new_right, StartVec[i].x_new_left);
		}

		M_big = StartVec[0].M;

		for (int i = 1; i < procs; i++) // найдем макс из остальных
		{
			if (M_big < StartVec[i].M)
			{
				M_big = StartVec[i].M;
			}
		}
		Lipschitz(M_big, m_small);// оценка параметра m

		 //подсчёт R
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			StartVec[i].R = R(m_small, StartVec[i].z_right, StartVec[i].z_left, StartVec[i].x_new_right, StartVec[i].x_new_left);
		}
		for (register int i = 0; i < procs; i++)//вставка результатов в общий список
		{
			Local_group.push_back(StartVec[i]);
		}
		Local_group.sort();//сортировка по возрастанию

						 //получение текущих границ
		curr_left = Local_group.back().x_new_left;
		curr_right = Local_group.back().x_new_right;

		localSteps++;

		while (true)//основной цикл
		{
			Pointer = --Local_group.end();
			for (int i = 0; i < procs * 2 - 2; i += 2) // хватаем proc_count R-ок (самых больших)
			{
				MainVec[i] = *Pointer;
				Local_group.erase(Pointer--);//остальные R-ки стираем (все группы)
			}
			MainVec[procs * 2 - 2] = *Pointer;
			Local_group.erase(Pointer);//.erase сотрет по указателю
			size_group = Local_group.size();

			// соответствующей  кучке больших R-ок находим соответствующие точки испытаний
#pragma omp parallel for shared(x_new)  num_threads(procs)
			for (int i = 0; i < procs; i++)
			{
				double coordinates_y[2];
				x_new[i] = New_x(MainVec[i * 2].x_new_right, MainVec[i * 2].x_new_left, MainVec[i * 2].z_right, MainVec[i * 2].z_left, m_small);

				evolve[i]->GetImage(x_new[i], coordinates_y);
				MainVec[i * 2 + 1].x_left = coordinates_y[0];
				MainVec[i * 2 + 1].y_left = coordinates_y[1];
				MainVec[i * 2 + 1].z_left = Func(coordinates_y, func);
			}

			for (int i = 0; i < procs * 2; i += 2)//сохраняем найденные точки
			{
				MainVec[i + 1].x_right = MainVec[i].x_right;
				MainVec[i + 1].y_right = MainVec[i].y_right;
				MainVec[i + 1].x_new_right = MainVec[i].x_new_right;
				MainVec[i + 1].z_right = MainVec[i].z_right;

				MainVec[i + 1].x_new_left = MainVec[i].x_new_right = x_new[i / 2];

				MainVec[i].x_right = MainVec[i + 1].x_left;
				MainVec[i].y_right = MainVec[i + 1].y_left;
				MainVec[i].z_right = MainVec[i + 1].z_left;
			}

			// считаем М
#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)
			{
				MainVec[i].M = M(MainVec[i].z_right, MainVec[i].z_left, MainVec[i].x_new_right, MainVec[i].x_new_left);
			}
			double M_max_array = MainVec[0].M;//найденная М одним потоком


	   //поиск наибольшего M, который не соответсвует текущей кучке(из M_max_all_M) (а он может быть большим, чем все M из кучки)
			for (int i = 1; i < procs * 2; i++)
			{
				if (M_max_array < MainVec[i].M)
				{
					M_max_array = MainVec[i].M;
				}
			}

			// поиск наибольшего M из кучки (из M_max_cogorta)
			if (size_group != 0)
			{
				Pointer = Local_group.begin();
				double tmp_M;
				double M_max_cogorta = (*Pointer++).M;

				while (Pointer != Local_group.end())
				{
					tmp_M = (*Pointer).M;
					if (M_max_cogorta < tmp_M) M_max_cogorta = tmp_M;

					Pointer++;
				}

				//поиск большего M среди найденных
				if (M_max_array < M_max_cogorta)
				{
					M_big = M_max_cogorta;
				}
				else
				{
					M_big = M_max_array;
				}

			}
			else
			{
				M_big = M_max_array;
			}


			//подсчёт m 
			Lipschitz(M_big, m_small);

#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)//считаем R для текущих отрезков
			{
				MainVec[i].R = R(m_small, MainVec[i].z_right, MainVec[i].z_left, MainVec[i].x_new_right, MainVec[i].x_new_left);
			}

			if (size_group != 0)
			{
				std::vector<GroupTwoDim> R_vec(Local_group.begin(), Local_group.end());
#pragma omp parallel for num_threads(procs)
				for (int i = 0; i < size_group; i++)//подсчёт R для отрезков  из очереди
				{
					R_vec[i].R = R(m_small, R_vec[i].z_right, R_vec[i].z_left, R_vec[i].x_new_right, R_vec[i].x_new_left);
				}
				int j = 0;
				for (std::list<GroupTwoDim>::iterator i = Local_group.begin(); i != Local_group.end(); i++)
				{
					i->R = R_vec[j++].R;
				}
				R_vec.clear();
			}
			for (int i = 0; i < procs * 2; i++)
			{
				Local_group.push_back(MainVec[i]);//кладем новые группы в общий список
			}
			//сбор данных 
			p.X.push_back(Local_group.back().x_right);
			p.Y.push_back(Local_group.back().y_right);

			Local_group.sort();
			localSteps++;
			if (localSteps > _Steps) { break; }//условие остановки
			curr_left = Local_group.back().x_new_left;
			curr_right = Local_group.back().x_new_right;
			if (abs(curr_left - curr_right) < eps) { break; }//условие остановки
		}
	}

	end = omp_get_wtime();
	//сбор данных 
	p.time = end - start;
	omp_set_dynamic(1);



	//сбор данных 
	p.X.push_back(Local_group.back().x_right);
	p.Y.push_back(Local_group.back().y_right);
	p.x = Local_group.back().x_right;
	p.y = Local_group.back().y_right;
	p.z = Local_group.back().z_right;

	// очистка
	StartVec.clear();
	MainVec.clear();
	Local_group.clear();
	evolve.clear();
	delete[]x_new;
	delete[]coordinates_y;


	return p;
}

bool operator<(const GroupTwoDim& one, const GroupTwoDim& two)
{
	if (one.R < two.R)
	{
		return true;
	}
	else
	{
		return false;
	}
};

bool operator>(const GroupTwoDim& one, const GroupTwoDim& two)
{
	if (one.R > two.R)
	{
		return true;
	}
	else
	{
		return false;
	}
};