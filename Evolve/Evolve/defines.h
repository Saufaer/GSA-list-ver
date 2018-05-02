#ifndef __DEFINES_H__
#define __DEFINES_H__


/* ======================================================================== *\
**  Constants                                                               **
\* ======================================================================== */
#define NOT_USE_CL
/// ������������ �����������
#define MAX_DIM 50
// ������������ ���������� ���������� ���������
#define MAX_NUM_MIN 20
// ���� ��������� ���������� �� �� GPU ����� �������������� double ����� float
#define CUDA_VALUE_DOUBLE_PRECISION
/// ���� ��������� ���������� �� �� �PU ����� �������������� double ����� float

#define CPU_VALUE_DOUBLE_PRECISION

/**
������ �������� ���������� �������, ������ �������� ��������� ������� ��������� �������.

�������� ������������� ������� � ���� ������ ��������� ����� ����������, ���� ������� ���������
���������� �������� �� ������ ������ (��������, ������� ����������).
*/
#define _M_ZERO_LEVEL 1e-12


/* ======================================================================== *\
**  Types                                                                   **
\* ======================================================================== */
/// ������������ ��� ������ ��� ���������� �� GPU
#ifdef CUDA_VALUE_DOUBLE_PRECISION
#define CUDA_VALUE double
#else
#define CUDA_VALUE float
#endif
/// ������������ ��� ������ ��� ���������� �� ����������� ����������(�� ��������),
#ifdef CPU_VALUE_DOUBLE_PRECISION
#define CPU_VALUE double
#else
#define CPU_VALUE float
#endif
/// ������ ����� �����
#define ints int*
/// ������ �������������� �����
#define doubles double*
/// ��� �����
#define FLAG bool

/* ======================================================================== *\
**  Defines                                                                 **
\* ======================================================================== */

/* ======================================================================== *\
**  ��������� ������������� � ������,  ��� ������ ������ ������������       **
\* ======================================================================== */

#ifdef _GPU_CUDA_ ///��� ���������� ��� CUDA

#define perm 1
#define SPECIFIER __shared__
#define ARRAY_SPECIFIER __constant__
#define concatenation cuda
#define F_DEVICE __device__
#define parameter_const
#define OBJECTIV_TYPE CUDA_VALUE
#define GET_FUNCTION_PARAMETERS OBJECTIV_TYPE* x, OBJECTIV_TYPE* f
#define FUNCTION_CALCOLATION_PREF (x)

#define GKLS_VARIABLES_SPECIFIER __shared__
//������������ �������� ���������� ��������
#ifdef CUDA_VALUE_DOUBLE_PRECISION
#define PRECISION(x) x##0
#else
#define PRECISION(x) x##f
#endif

#else ///��� ���������� ��� VS

#define perm 0
#define SPECIFIER extern
#define ARRAY_SPECIFIER extern
#define concatenation
#define F_DEVICE inline
#define parameter_const const
#define OBJECTIV_TYPE CPU_VALUE
#define GET_FUNCTION_PARAMETERS tFunction* f
#define FUNCTION_CALCOLATION_PREF

#define GKLS_VARIABLES_SPECIFIER extern
//������������ �������� ���������� ��������
#ifdef CPU_VALUE_DOUBLE_PRECISION
#define PRECISION(x) x##0
#else
#define PRECISION(x) x##f
#endif

#endif

/* ======================================================================== *\
**  ���������                                                               **
\* ======================================================================== */
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? a : b)
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? a : b)
#endif

/// ���������� ��� �����
#define CAT(x, y) x##y
/// ���������� ������ �����
#define CAT4(a, b, c, d) a##b##c##d
/// ���������� ��� ����� � ���������������
#define CONCATENATION2(name, console) CAT(console,name)
/** ��������� ������� concatenation
(��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
*/
#define CONCATENATION(name) CONCATENATION2(name, concatenation)
/// ���������� ����� � �����������
#define CAT_COM(x) //##x

/// ��������� ������� P � ����� ����, ��� ������������� ������������ EParameterType
#define ParType(type) P##type
/// ��������� ������� link � �����
#define LinkParameter(name) link##name
/// ��������� ������� com � �����
#define ComParameter(name) com##name
/// ��������� ������� help � �����
#define HelpParameter(name) help##name
/// ��������� ������� inc � �����
#define IncParameter(name) inc##name
/// ��������� ����� bar � ������ ���� char*
#define make_str(bar) # bar
/// ��������� ������� IS_ � �����
#define IsChange(name) IS_##name

/// �������������� �������� �� ������ TParameters
#define InitParameter(type, name, defVal, com, help, sizeVal)                            \
IsChange(name) = false;                                                                  \
Inc(ParType(name), ParType(type), LinkParameter(name), com,                              \
  HelpParameter(name), help, (void*)(&name), make_str(defVal),                           \
  make_str(name), &IsChange(name), IncParameter(name), make_str(name), sizeVal);

#define OWNER_NAME Owner

/// �������������� �������� �� ������ ����������
#define InitOption(name, defVal, com, help, sizeVal)                           \
  InitializationOption((TBaseProperty<OWNER_NAME>*)(&name), make_str(name), make_str(defVal), com, help, sizeVal);

/* ======================================================================== *\
**  �������� ��� ���������� ����������                                      **
\* ======================================================================== */

/** ������� ���������� ����������
name - ��� ����������
type - ��� ����������
specifier - ������������
*/
#define VARIABLES(name, type, specifier) specifier type name
/** ���������� ���������� �� �������������� GKLS_VARIABLES_SPECIFIER
(��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
*/
#define GKLS_VARIABLES(name, type) VARIABLES(name, type, GKLS_VARIABLES_SPECIFIER)
/** ���������� ���������� �� �������������� SPECIFIER
(��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
*/
#define CONSTANT_VARIABLES(name, type) SPECIFIER type name
/// ��������� ������� ��� ������������� �������� �� ���������
#define NEW_FUNC_DEF(name) CONCATENATION(name##_func_def())
/// ������� ����������� ������
#define STATIC_ARRAY(type, name, count) type name[count]
/** ������� ����������� ������ �� �������������� ARRAY_SPECIFIER
(��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
*/
#define NEW_ARRAY(type, name, count) ARRAY_SPECIFIER STATIC_ARRAY(type, name, count)
/** ������� ����������� ������ �������� MAX_DIM ���� OBJECTIV_TYPE �� �������������� ARRAY_SPECIFIER
(��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
*/
#define NEW_ARRAY_MAX_SIZE(name) NEW_ARRAY(OBJECTIV_TYPE, name, MAX_DIM)
/** ������� ���������� ���� OBJECTIV_TYPE �� �������������� SPECIFIER
(��.��������� ������������� � ������, ��� ������ ������ ������������)
*/
#define FLOAT_VARIABLES(name) SPECIFIER OBJECTIV_TYPE name


/**���������� ���������� ����������� � ����� � �������, ����� ���� ��������������� � ����������,
!!!� problems.cpp ��������� ��������� ����������: type name;!!!
type - ���
name - ���
val - �������� �� ���������
*/
#define NEW_VARIABLES(type, name, val) CONSTANT_VARIABLES(name, type);\
 F_DEVICE type NEW_FUNC_DEF(name)\
  {\
  name = val;\
return name;\
}\

/// ��������� ���������� ��� �������� ��������� ������ TParameters
#define CreateParameter(type, name)       \
public: type name;                        \
public: FLAG IsChange(name);              \
protected: int IncParameter(name);        \
protected: EParameterType ParType(name);  \
protected: string LinkParameter(name);    \
protected: string HelpParameter(name);

/// ������� ��������������� ��� ������� ����� ������
#define BasicMethods(ClassType, Type)                                             \
  virtual void operator =(Type data) {TTypedProperty<Type, Owner>::operator=(data);}                \
  virtual void operator=(ClassType<Owner>& data) {TParameterProperty<Type, Owner>::operator=(data);}     \
  virtual void operator = (ClassType<Owner>* data) {TParameterProperty<Type, Owner>::operator=(data);}  \
  virtual string ToString() {return operator string();}                           \
  virtual void FromString(string val) {operator = (val);}                           \
  virtual void operator = (TBaseProperty<Owner>& data) { TParameterProperty<Type, Owner>::operator=((TParameterProperty<Type, Owner>*)(& data)); } \
  virtual void operator=(char* data) {operator = (string(data));}

/* ======================================================================== *\
**  �������� ��� ���������� �������                                         **
\* ======================================================================== */

/** ������� ������� � ������� concatenation##name, ������������ �������� ���� type
(�������� concatenation ��. ��������� ������������� � ������,  ��� ������ ������ ������������  )
�� �������������� console,
���������� parameter ���� ������ type � ������������� parameterconst
*/
#define OBJECTIV_FUNCTION2(name, type, parameter, parameterconst, console) console type CONCATENATION(name)(parameterconst type parameter[])
/** ������� ������� � ������� concatenation##name, ������������ �������� ���� type
�� �������������� F_DEVICE,
���������� parameter ���� ������ type � ������������� parameter_const
(��.��������� ������������� � ������, ��� ������ ������ ������������)
*/
#define OBJECTIV_FUNCTION(name, type, parameter) OBJECTIV_FUNCTION2(name, type, parameter, parameter_const, F_DEVICE)
/** ������� ������� � ������� concatenation##name, ������������ �������� ���� type
�� �������������� console,
����� ����������� parameter1, parameter2 ���� ������ type � ������������� parameterconst
*/
#define P2_FUNCTION2(name, type, parameter1, parameter2, parameterconst, console) console type CONCATENATION(name)(parameterconst type* parameter1, parameterconst type* parameter2)
/** ������� ������� � ������� concatenation##namen, ������������ �������� ���� type
�� �������������� F_DEVICE,
����� ����������� parameter1, parameter2 ���� ������ type � ������������� parameter_const
*/
#define P2_FUNCTION(name, type, parameter1, parameter2) P2_FUNCTION2(name, type, parameter1, parameter2, parameter_const, F_DEVICE)
//////////////////////////////////////////////////////////////////////////

/** �������  � ������� concatenation##name
���������� �������� ���� OBJECTIV_TYPE, �� �������������� F_DEVICE,
parameter - ��� ���������, ���� OBJECTIV_TYPE*  ������������� parameterconst
(��.��������� ������������� � ������, ��� ������ ������ ������������)
*/
#define HYBRID_OBJECTIV_FUNCTION(name, parameter) OBJECTIV_FUNCTION(name, OBJECTIV_TYPE, parameter)
/** �������  � ������� concatenation##objfn
���������� �������� ���� OBJECTIV_TYPE, �� �������������� F_DEVICE,
� ���������� x , ���� OBJECTIV_TYPE*  ������������� parameterconst
(��.��������� ������������� � ������, ��� ������ ������ ������������)
*/
#define OBJECTIV HYBRID_OBJECTIV_FUNCTION(objfn, x)
/** ������� �����������, ���� ��� concatenation##constraint##num,
���������� �������� ���� OBJECTIV_TYPE � �� �������������� F_DEVICE,
� ���������� x , ���� OBJECTIV_TYPE*  ������������� parameterconst
(��.��������� ������������� � ������, ��� ������ ������ ������������)
num - ����� �������
*/
#define CONSTRAINT(num) HYBRID_OBJECTIV_FUNCTION(CAT(constraint,num), x)

/** ������� �� ��� ���������� ���������� �������� ���� OBJECTIV_TYPE,
name - ���
parameter1 � parameter2 - ����� ����������, ���� OBJECTIV_TYPE*
���������� �������� ���� OBJECTIV_TYPE
*/
#define HYBRID_P2_FUNCTION2(name, parameter1, parameter2) P2_FUNCTION(name, OBJECTIV_TYPE, parameter1, parameter2)
/** �������� ������� ��� ������������� � CUDA,
���������:
type - ��� ������������� ��������
name - ��� �������
*/
#define NEW_FUNCTION(type, name) F_DEVICE type CONCATENATION(name)
/// ������ ������ ������� ������
#define GET_FUNCTION NEW_FUNCTION(void, GetFunction)(GET_FUNCTION_PARAMETERS)
/// ������� �������
#define OBJ CONCATENATION(objfn) FUNCTION_CALCOLATION_PREF
/// ������� �����������
#define CON(num) CONCATENATION(CAT(constraint,num)) FUNCTION_CALCOLATION_PREF

/**���������� Set � Get ������� c ������
Get##N ���������� ��� T
� Set##N �������� value ���� T
*/
#define PROPERTY(T, N)     \
  T Get ## N() const;     \
  void Set ## N(T value);

/* ======================================================================== *\
**  ������ ���������                                                        **
\* ======================================================================== */

#define PI PRECISION(3.14159265359)







#endif

