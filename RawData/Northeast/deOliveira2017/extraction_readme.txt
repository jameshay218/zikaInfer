Data included here is, all for Northeast Brazil:
1. Digitised confirmed and predicted microcephaly incidence per 10000 births
2. Digitised suspected ZIKV per 10000 pregnant women
3. Correct number of infected pregnant women
4. Correct total number of microcephaly cases

Given that 2. represents incidence per 10000 women, and we have actual incidence, we can calculate the number of pregnant women that were assumed in each month as:
number_cases/inc_per_person

This can then be used to infer the number of monthly births, as we are told that number of pregnant women is calculated from number of reported live births as:

pregnant_N = (births*9) +(births*0.2)*1.5
	   = births*9.3
births=pregnant_N/9.3

We can then infer the number of microcephaly cases per month as:
microceph_cases = births*microceph_inc_per_birth

I then adjusted these estimated monthly cases by scaling to give 1487 total cases

**Note** that the total population size for the Northeast, N, is still the same, as the reporting proportion parameter will mop up the fact that only a certain proportion of the total population are pregnant women. We have to interpret this parameter as "probability that an individual was a pregnant woman AND was reported to the healthcare authorities, GIVEN that they were infected"
