﻿#pragma once

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#include "Search.h"
#include <omp.h>
namespace F {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Сводка для FormAGP
	/// </summary>
	public ref class FormAGP : public System::Windows::Forms::Form
	{
	public:
		FormAGP(void)
		{
			InitializeComponent();
			//
			//TODO: добавьте код конструктора
			//
		}

	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~FormAGP()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::TextBox^  BoxminZ;
	private: System::Windows::Forms::Button^  Butcalc;
	protected:

	protected:


	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::TextBox^  BoxmaxN;
	private: System::Windows::Forms::TextBox^  Boxeps;




	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::TextBox^  BoxminX;


	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  Boxsteps;

	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::TextBox^  Boxtime;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::TextBox^  threadsBox;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^  chart1;



	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::TextBox^  textBox3;
	private: System::Windows::Forms::TextBox^  textBox4;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::Label^  label10;
	private: System::Windows::Forms::TextBox^  textBox5;
	private: System::Windows::Forms::Label^  label11;
	private: System::Windows::Forms::Label^  label12;
	private: System::Windows::Forms::Label^  label13;
	private: System::Windows::Forms::TextBox^  textBox6;
	private: System::Windows::Forms::TextBox^  textBox7;
	private: System::Windows::Forms::Label^  label14;
	private: System::Windows::Forms::RadioButton^  radioButton1;
	private: System::Windows::Forms::RadioButton^  radioButton2;
	private: System::Windows::Forms::TabControl^  tabControl1;
	private: System::Windows::Forms::TabPage^  tabPage1;
	private: System::Windows::Forms::TabPage^  tabPage2;
	private: System::Windows::Forms::DataGridView^  dataGridView1;



	private: System::Windows::Forms::GroupBox^  groupBox1;
	private: System::Windows::Forms::GroupBox^  groupBox2;
	private: System::Windows::Forms::GroupBox^  groupBox3;
	private: System::Windows::Forms::TabPage^  tabPage3;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  Num;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  X;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  Z;











	protected:

	private:
		/// <summary>
		/// Обязательная переменная конструктора.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Требуемый метод для поддержки конструктора — не изменяйте 
		/// содержимое этого метода с помощью редактора кода.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataVisualization::Charting::ChartArea^  chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^  legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^  series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^  series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^  series3 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			this->BoxminZ = (gcnew System::Windows::Forms::TextBox());
			this->Butcalc = (gcnew System::Windows::Forms::Button());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->BoxmaxN = (gcnew System::Windows::Forms::TextBox());
			this->Boxeps = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->BoxminX = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->Boxsteps = (gcnew System::Windows::Forms::TextBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->Boxtime = (gcnew System::Windows::Forms::TextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->threadsBox = (gcnew System::Windows::Forms::TextBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->textBox6 = (gcnew System::Windows::Forms::TextBox());
			this->textBox7 = (gcnew System::Windows::Forms::TextBox());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->tabControl1 = (gcnew System::Windows::Forms::TabControl());
			this->tabPage3 = (gcnew System::Windows::Forms::TabPage());
			this->groupBox1 = (gcnew System::Windows::Forms::GroupBox());
			this->tabPage1 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage2 = (gcnew System::Windows::Forms::TabPage());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->Num = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->X = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Z = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->groupBox2 = (gcnew System::Windows::Forms::GroupBox());
			this->groupBox3 = (gcnew System::Windows::Forms::GroupBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->tabControl1->SuspendLayout();
			this->tabPage3->SuspendLayout();
			this->groupBox1->SuspendLayout();
			this->tabPage1->SuspendLayout();
			this->tabPage2->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->groupBox2->SuspendLayout();
			this->groupBox3->SuspendLayout();
			this->SuspendLayout();
			// 
			// BoxminZ
			// 
			this->BoxminZ->Enabled = false;
			this->BoxminZ->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->BoxminZ->Location = System::Drawing::Point(10, 49);
			this->BoxminZ->Name = L"BoxminZ";
			this->BoxminZ->Size = System::Drawing::Size(215, 29);
			this->BoxminZ->TabIndex = 0;
			this->BoxminZ->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// Butcalc
			// 
			this->Butcalc->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Butcalc->Location = System::Drawing::Point(7, 411);
			this->Butcalc->Name = L"Butcalc";
			this->Butcalc->Size = System::Drawing::Size(303, 49);
			this->Butcalc->TabIndex = 1;
			this->Butcalc->Text = L"Расчёт";
			this->Butcalc->UseVisualStyleBackColor = true;
			this->Butcalc->Click += gcnew System::EventHandler(this, &FormAGP::Butcalc_Click);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(6, 22);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(199, 24);
			this->label1->TabIndex = 2;
			this->label1->Text = L"Найденный минимум";
			// 
			// BoxmaxN
			// 
			this->BoxmaxN->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->BoxmaxN->Location = System::Drawing::Point(13, 160);
			this->BoxmaxN->Name = L"BoxmaxN";
			this->BoxmaxN->Size = System::Drawing::Size(215, 29);
			this->BoxmaxN->TabIndex = 3;
			this->BoxmaxN->Text = L"800";
			this->BoxmaxN->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// Boxeps
			// 
			this->Boxeps->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Boxeps->Location = System::Drawing::Point(13, 230);
			this->Boxeps->Name = L"Boxeps";
			this->Boxeps->Size = System::Drawing::Size(215, 29);
			this->Boxeps->TabIndex = 4;
			this->Boxeps->Text = L"0,0001";
			this->Boxeps->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(12, 133);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(291, 24);
			this->label2->TabIndex = 5;
			this->label2->Text = L"Максимальное число итераций";
			// 
			// BoxminX
			// 
			this->BoxminX->Enabled = false;
			this->BoxminX->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->BoxminX->Location = System::Drawing::Point(11, 108);
			this->BoxminX->Name = L"BoxminX";
			this->BoxminX->Size = System::Drawing::Size(215, 29);
			this->BoxminX->TabIndex = 6;
			this->BoxminX->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(6, 81);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(220, 24);
			this->label3->TabIndex = 7;
			this->label3->Text = L"Найденная координата";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label4->Location = System::Drawing::Point(11, 203);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(244, 24);
			this->label4->TabIndex = 8;
			this->label4->Text = L"Допустимая погрешность";
			// 
			// Boxsteps
			// 
			this->Boxsteps->Enabled = false;
			this->Boxsteps->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Boxsteps->Location = System::Drawing::Point(11, 176);
			this->Boxsteps->Name = L"Boxsteps";
			this->Boxsteps->Size = System::Drawing::Size(215, 29);
			this->Boxsteps->TabIndex = 9;
			this->Boxsteps->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label5->Location = System::Drawing::Point(7, 149);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(282, 24);
			this->label5->TabIndex = 10;
			this->label5->Text = L"Число проведенных итераций";
			// 
			// Boxtime
			// 
			this->Boxtime->Enabled = false;
			this->Boxtime->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->Boxtime->Location = System::Drawing::Point(11, 235);
			this->Boxtime->Name = L"Boxtime";
			this->Boxtime->Size = System::Drawing::Size(215, 29);
			this->Boxtime->TabIndex = 11;
			this->Boxtime->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label6->Location = System::Drawing::Point(7, 208);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(269, 24);
			this->label6->TabIndex = 12;
			this->label6->Text = L"Затраченное время (секунд)";
			// 
			// threadsBox
			// 
			this->threadsBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->threadsBox->Location = System::Drawing::Point(11, 51);
			this->threadsBox->Name = L"threadsBox";
			this->threadsBox->Size = System::Drawing::Size(71, 29);
			this->threadsBox->TabIndex = 13;
			this->threadsBox->Text = L"1";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label7->Location = System::Drawing::Point(7, 24);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(75, 24);
			this->label7->TabIndex = 14;
			this->label7->Text = L"Потоки";
			// 
			// chart1
			// 
			this->chart1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			chartArea1->BackImageAlignment = System::Windows::Forms::DataVisualization::Charting::ChartImageAlignmentStyle::Top;
			chartArea1->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea1);
			legend1->Name = L"Legend1";
			this->chart1->Legends->Add(legend1);
			this->chart1->Location = System::Drawing::Point(3, 6);
			this->chart1->Name = L"chart1";
			series1->BorderWidth = 2;
			series1->ChartArea = L"ChartArea1";
			series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Spline;
			series1->Color = System::Drawing::Color::DarkTurquoise;
			series1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			series1->Legend = L"Legend1";
			series1->Name = L"Track";
			series2->BorderWidth = 2;
			series2->ChartArea = L"ChartArea1";
			series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::FastPoint;
			series2->Color = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(128)), static_cast<System::Int32>(static_cast<System::Byte>(128)),
				static_cast<System::Int32>(static_cast<System::Byte>(255)));
			series2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			series2->Legend = L"Legend1";
			series2->Name = L"Points";
			series2->YValuesPerPoint = 2;
			series3->BorderWidth = 4;
			series3->ChartArea = L"ChartArea1";
			series3->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::FastPoint;
			series3->Color = System::Drawing::Color::Black;
			series3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			series3->LabelBorderWidth = 4;
			series3->Legend = L"Legend1";
			series3->MarkerColor = System::Drawing::Color::Red;
			series3->MarkerSize = 10;
			series3->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Cross;
			series3->Name = L"Result";
			series3->YValuesPerPoint = 4;
			this->chart1->Series->Add(series1);
			this->chart1->Series->Add(series2);
			this->chart1->Series->Add(series3);
			this->chart1->Size = System::Drawing::Size(769, 413);
			this->chart1->TabIndex = 15;
			this->chart1->Text = L"chart1";
			// 
			// textBox1
			// 
			this->textBox1->BackColor = System::Drawing::SystemColors::Menu;
			this->textBox1->Enabled = false;
			this->textBox1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox1->Location = System::Drawing::Point(12, 6);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(30, 29);
			this->textBox1->TabIndex = 19;
			this->textBox1->Text = L"2";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label8->Location = System::Drawing::Point(48, 9);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(37, 24);
			this->label8->TabIndex = 20;
			this->label8->Text = L"Sin";
			// 
			// textBox2
			// 
			this->textBox2->BackColor = System::Drawing::SystemColors::Menu;
			this->textBox2->Enabled = false;
			this->textBox2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox2->Location = System::Drawing::Point(91, 6);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(35, 29);
			this->textBox2->TabIndex = 21;
			this->textBox2->Text = L"10";
			// 
			// textBox3
			// 
			this->textBox3->BackColor = System::Drawing::SystemColors::Menu;
			this->textBox3->Enabled = false;
			this->textBox3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox3->Location = System::Drawing::Point(189, 6);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(28, 29);
			this->textBox3->TabIndex = 22;
			this->textBox3->Text = L"3";
			// 
			// textBox4
			// 
			this->textBox4->BackColor = System::Drawing::SystemColors::Menu;
			this->textBox4->Enabled = false;
			this->textBox4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox4->Location = System::Drawing::Point(272, 6);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(38, 29);
			this->textBox4->TabIndex = 23;
			this->textBox4->Text = L"11";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label9->Location = System::Drawing::Point(223, 9);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(43, 24);
			this->label9->TabIndex = 24;
			this->label9->Text = L"Cos";
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label10->Location = System::Drawing::Point(162, 11);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(21, 24);
			this->label10->TabIndex = 25;
			this->label10->Text = L"+";
			// 
			// textBox5
			// 
			this->textBox5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox5->Location = System::Drawing::Point(13, 96);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(101, 29);
			this->textBox5->TabIndex = 26;
			this->textBox5->Text = L"2";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label11->Location = System::Drawing::Point(11, 69);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(110, 24);
			this->label11->TabIndex = 27;
			this->label11->Text = L"Параметр r";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label12->Location = System::Drawing::Point(132, 9);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(24, 24);
			this->label12->TabIndex = 28;
			this->label12->Text = L"X";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label13->Location = System::Drawing::Point(316, 9);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(24, 24);
			this->label13->TabIndex = 29;
			this->label13->Text = L"X";
			// 
			// textBox6
			// 
			this->textBox6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox6->Location = System::Drawing::Point(55, 33);
			this->textBox6->Name = L"textBox6";
			this->textBox6->Size = System::Drawing::Size(44, 29);
			this->textBox6->TabIndex = 30;
			this->textBox6->Text = L"2";
			// 
			// textBox7
			// 
			this->textBox7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox7->Location = System::Drawing::Point(105, 33);
			this->textBox7->Name = L"textBox7";
			this->textBox7->Size = System::Drawing::Size(42, 29);
			this->textBox7->TabIndex = 31;
			this->textBox7->Text = L"8";
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label14->Location = System::Drawing::Point(9, 36);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(40, 24);
			this->label14->TabIndex = 32;
			this->label14->Text = L"X ∈";
			// 
			// radioButton1
			// 
			this->radioButton1->AutoSize = true;
			this->radioButton1->Checked = true;
			this->radioButton1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->radioButton1->Location = System::Drawing::Point(111, 25);
			this->radioButton1->Name = L"radioButton1";
			this->radioButton1->Size = System::Drawing::Size(174, 28);
			this->radioButton1->TabIndex = 33;
			this->radioButton1->TabStop = true;
			this->radioButton1->Text = L"Характеристики";
			this->radioButton1->UseVisualStyleBackColor = true;
			// 
			// radioButton2
			// 
			this->radioButton2->AutoSize = true;
			this->radioButton2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->radioButton2->Location = System::Drawing::Point(111, 57);
			this->radioButton2->Name = L"radioButton2";
			this->radioButton2->Size = System::Drawing::Size(187, 28);
			this->radioButton2->TabIndex = 34;
			this->radioButton2->Text = L"Простое деление";
			this->radioButton2->UseVisualStyleBackColor = true;
			// 
			// tabControl1
			// 
			this->tabControl1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->tabControl1->Controls->Add(this->tabPage3);
			this->tabControl1->Controls->Add(this->tabPage1);
			this->tabControl1->Controls->Add(this->tabPage2);
			this->tabControl1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->tabControl1->Location = System::Drawing::Point(357, 9);
			this->tabControl1->Name = L"tabControl1";
			this->tabControl1->SelectedIndex = 0;
			this->tabControl1->Size = System::Drawing::Size(786, 451);
			this->tabControl1->TabIndex = 35;
			// 
			// tabPage3
			// 
			this->tabPage3->Controls->Add(this->groupBox1);
			this->tabPage3->Location = System::Drawing::Point(4, 27);
			this->tabPage3->Name = L"tabPage3";
			this->tabPage3->Size = System::Drawing::Size(778, 420);
			this->tabPage3->TabIndex = 2;
			this->tabPage3->Text = L"Результаты";
			this->tabPage3->UseVisualStyleBackColor = true;
			// 
			// groupBox1
			// 
			this->groupBox1->Controls->Add(this->label1);
			this->groupBox1->Controls->Add(this->BoxminZ);
			this->groupBox1->Controls->Add(this->label3);
			this->groupBox1->Controls->Add(this->BoxminX);
			this->groupBox1->Controls->Add(this->label5);
			this->groupBox1->Controls->Add(this->Boxsteps);
			this->groupBox1->Controls->Add(this->label6);
			this->groupBox1->Controls->Add(this->Boxtime);
			this->groupBox1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->groupBox1->Location = System::Drawing::Point(18, 16);
			this->groupBox1->Name = L"groupBox1";
			this->groupBox1->Size = System::Drawing::Size(720, 311);
			this->groupBox1->TabIndex = 36;
			this->groupBox1->TabStop = false;
			// 
			// tabPage1
			// 
			this->tabPage1->Controls->Add(this->chart1);
			this->tabPage1->Location = System::Drawing::Point(4, 27);
			this->tabPage1->Name = L"tabPage1";
			this->tabPage1->Padding = System::Windows::Forms::Padding(3);
			this->tabPage1->Size = System::Drawing::Size(778, 420);
			this->tabPage1->TabIndex = 0;
			this->tabPage1->Text = L"График";
			this->tabPage1->UseVisualStyleBackColor = true;
			// 
			// tabPage2
			// 
			this->tabPage2->Controls->Add(this->dataGridView1);
			this->tabPage2->Location = System::Drawing::Point(4, 27);
			this->tabPage2->Name = L"tabPage2";
			this->tabPage2->Padding = System::Windows::Forms::Padding(3);
			this->tabPage2->Size = System::Drawing::Size(778, 420);
			this->tabPage2->TabIndex = 1;
			this->tabPage2->Text = L"Таблица точек поиска";
			this->tabPage2->UseVisualStyleBackColor = true;
			// 
			// dataGridView1
			// 
			this->dataGridView1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(3) {
				this->Num, this->X,
					this->Z
			});
			this->dataGridView1->Location = System::Drawing::Point(-4, 0);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->Size = System::Drawing::Size(782, 417);
			this->dataGridView1->TabIndex = 0;
			// 
			// Num
			// 
			this->Num->HeaderText = L"Number";
			this->Num->Name = L"Num";
			// 
			// X
			// 
			this->X->HeaderText = L"X";
			this->X->Name = L"X";
			this->X->Width = 200;
			// 
			// Z
			// 
			this->Z->HeaderText = L"Z";
			this->Z->Name = L"Z";
			this->Z->Width = 200;
			// 
			// groupBox2
			// 
			this->groupBox2->Controls->Add(this->label2);
			this->groupBox2->Controls->Add(this->BoxmaxN);
			this->groupBox2->Controls->Add(this->Boxeps);
			this->groupBox2->Controls->Add(this->label4);
			this->groupBox2->Controls->Add(this->label14);
			this->groupBox2->Controls->Add(this->textBox7);
			this->groupBox2->Controls->Add(this->label11);
			this->groupBox2->Controls->Add(this->textBox6);
			this->groupBox2->Controls->Add(this->textBox5);
			this->groupBox2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->groupBox2->Location = System::Drawing::Point(7, 42);
			this->groupBox2->Name = L"groupBox2";
			this->groupBox2->Size = System::Drawing::Size(306, 267);
			this->groupBox2->TabIndex = 37;
			this->groupBox2->TabStop = false;
			this->groupBox2->Text = L"Входные параметры";
			// 
			// groupBox3
			// 
			this->groupBox3->Controls->Add(this->label7);
			this->groupBox3->Controls->Add(this->threadsBox);
			this->groupBox3->Controls->Add(this->radioButton1);
			this->groupBox3->Controls->Add(this->radioButton2);
			this->groupBox3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->groupBox3->Location = System::Drawing::Point(4, 315);
			this->groupBox3->Name = L"groupBox3";
			this->groupBox3->Size = System::Drawing::Size(306, 90);
			this->groupBox3->TabIndex = 38;
			this->groupBox3->TabStop = false;
			this->groupBox3->Text = L"Способ распареллеливания";
			// 
			// FormAGP
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1155, 502);
			this->Controls->Add(this->groupBox3);
			this->Controls->Add(this->groupBox2);
			this->Controls->Add(this->tabControl1);
			this->Controls->Add(this->label13);
			this->Controls->Add(this->label12);
			this->Controls->Add(this->label10);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->Butcalc);
			this->Name = L"FormAGP";
			this->Text = L"AGP";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->tabControl1->ResumeLayout(false);
			this->tabPage3->ResumeLayout(false);
			this->groupBox1->ResumeLayout(false);
			this->groupBox1->PerformLayout();
			this->tabPage1->ResumeLayout(false);
			this->tabPage2->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->groupBox2->ResumeLayout(false);
			this->groupBox2->PerformLayout();
			this->groupBox3->ResumeLayout(false);
			this->groupBox3->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

	private: System::Void Butcalc_Click(System::Object^  sender, System::EventArgs^  e) {



		if (radioButton1->Checked == 1)
		{

			Pointer p1;
			double left = Convert::ToDouble(textBox6->Text);
			double right = Convert::ToDouble(textBox7->Text);
			double r = Convert::ToDouble(textBox5->Text);
			double eps = Convert::ToDouble(Boxeps->Text);

			int n_max = Convert::ToInt32(BoxmaxN->Text);
			int threads = Convert::ToInt32(threadsBox->Text);

			Search ptest1(left, right, r, threads);

			p1 = ptest1.Ch_SearchMin(eps, n_max, threads);

			chart1->ChartAreas["ChartArea1"]->AxisX->Minimum = left - 0.5;
			chart1->ChartAreas["ChartArea1"]->AxisX->Maximum = right + 0.5;



			chart1->ChartAreas["ChartArea1"]->AxisX->Title = "X";
			chart1->ChartAreas["ChartArea1"]->AxisY->Title = "F(X)";

			chart1->Series["Track"]->Points->Clear();
			chart1->Series["Result"]->Points->Clear();

			double minOY = ptest1.Func(0);
			for (double i = left; i <= right; i += 0.1)
			{
				double Z = ptest1.Func(i);
				chart1->Series["Track"]->Points->AddXY(i, Z);
				if (Z < minOY)
				{
					minOY = Z;
				}
			}
			chart1->Series["Result"]->Points->AddXY(p1.x, p1.z);

			chart1->Series["Points"]->Points->Clear();
			dataGridView1->Rows->Clear();
			for (int i = 0; i < p1.X.size(); i++)
			{

				chart1->Series["Points"]->Points->AddXY(p1.X[i], minOY - 1);
				dataGridView1->Rows->Insert(i, i, p1.X[i].ToString(), ptest1.Func(p1.X[i]).ToString());
			}

			Boxsteps->Text = Convert::ToString(p1.steps*threads);
			BoxminZ->Text = Convert::ToString(p1.z);
			BoxminX->Text = Convert::ToString(p1.x);
			Boxtime->Text = Convert::ToString(p1.time);



		}


		if (radioButton2->Checked == 1)
		{


			Pointer p;
			double left = Convert::ToDouble(textBox6->Text);
			double right = Convert::ToDouble(textBox7->Text);
			double r = Convert::ToDouble(textBox5->Text);
			double eps = Convert::ToDouble(Boxeps->Text);

			int n_max = Convert::ToInt32(BoxmaxN->Text);
			int threads = Convert::ToInt32(threadsBox->Text);

			Search ptest(left, right, r, threads);

			p = ptest.Simple_Par_searchMin(left, right, n_max, eps, threads);

			chart1->ChartAreas["ChartArea1"]->AxisX->Minimum = left - 0.5;
			chart1->ChartAreas["ChartArea1"]->AxisX->Maximum = right + 0.5;

			chart1->Series["Track"]->Points->Clear();
			chart1->Series["Result"]->Points->Clear();
			dataGridView1->Rows->Clear();
			double minOY = ptest.Func(0);
			for (double i = left; i <= right; i += 0.1)
			{
				double Z = ptest.Func(i);
				chart1->Series["Track"]->Points->AddXY(i, Z);
				if (Z < minOY)
				{
					minOY = Z;
				}
			}
			chart1->Series["Result"]->Points->AddXY(p.x, p.z);

			chart1->Series["Points"]->Points->Clear();

			for (int i = 0; i < p.X.size(); i++)
			{

				chart1->Series["Points"]->Points->AddXY(p.X[i], minOY - 1);
				dataGridView1->Rows->Insert(i, i, p.X[i].ToString(), ptest.Func(p.X[i]).ToString());

			}

			Boxsteps->Text = Convert::ToString(p.steps);
			BoxminZ->Text = Convert::ToString(p.z);
			BoxminX->Text = Convert::ToString(p.x);
			Boxtime->Text = Convert::ToString(p.time);

		}
	};
	};
}
