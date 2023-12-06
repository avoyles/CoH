// interface to C-code, MT19937

/**************************************/
/*      mt19937ar.c                   */
/**************************************/
#ifdef __cplusplus
extern "C" {
#endif
void init_by_array(unsigned long [], unsigned long);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
#ifdef __cplusplus
}
#endif
