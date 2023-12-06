// formatting of output values

/*** lower limit of output values, enforce zero */
static const double output_eps = 1.0e-99;


static string cline="#          ";
static string blank="           ";
static string dashl=" ----------";
static const int DisplayWidth = 99;

void outSectionHead (const char *);
double outRetrieveLabE (void);

static inline double lowfilter(double x)
{ if(fabs(x) < output_eps) return 0.0; else return(x); }

static inline void outVal(double x)
{ cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4) << setw(11) << x; }

static inline void outVal(int w, double x)
{ cout.setf(ios::scientific, ios::floatfield);
  int p = w - 8;
  if(p >= 1){ cout << setprecision(w-8) << setw(w) << x; }
  else      { cout << setw(w) << x; }}

static inline void outVal(int w, int x)
{ cout << dec << setw(w) << x; }

static inline void outVal(int w, unsigned int x)
{ cout << dec << setw(w) << x; }

static inline void outVal(int w, int p, double x)
{ cout.setf(ios::fixed, ios::floatfield);
  cout << setw(w) << setprecision(p) << x; }

static inline void nl()
{  cout << endl; }

