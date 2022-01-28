noinst_LTLIBRARIES = libBMDS.la 

AM_CPPFLAGS= -I$(top_srcdir)/code_base -I$(top_srcdir)/include  -I$(EIGEN_INCLUDE)
AM_CPPFLAGS+=  -O3 -funroll-loops -std=c++11
AM_CPPFLAGS+= -fopenmp -Wno-ignored-attributes -DEIGEN_NO_DEBUG -march=native

libBMDS_la_SOURCES=  DichLogProbitBMD_NC.cpp lognormal_HILL_NC.cpp \
		        bmdStruct.cpp dichotomous_entry_code.cpp  lognormalModels.cpp \
			DichProbitBMD_NC.cpp lognormal_POLYNOMIAL_NC.cpp \
			DichQlinearBMD_NC.cpp  lognormal_POWER_NC.cpp \
			DichWeibullBMD_NC.cpp  normal_EXP_NC.cpp \
			continuous_clean_aux.cpp   normal_HILL_NC.cpp normal_FUNL_NC.cpp\
			continuous_entry_code.cpp gradient.cpp  normalModels.cpp \
			continuous_model_functions.cpp  helperfunctions.cpp normal_POLYNOMIAL_NC.cpp \
			DichGammaBMD_NC.cpp IDPrior.cpp  normal_POWER_NC.cpp \
			DichLogisticBMD_NC.cpp lognormal_EXP_NC.cpp  normal_likelihoods.cpp stdafx.cpp





