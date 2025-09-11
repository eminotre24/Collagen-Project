from scipy import odr

def lin_mdl(beta, data):
    m, c = beta
    return m * data + c

def get_odr(x, y, beta0):
    fun = odr.Model(lin_mdl)
    data = odr.Data(x, y)
    odr_run = odr.ODR(data, fun, beta0=beta0) # Guess Parameters
    res = odr_run.run()
    m_res, b_res = res.beta
    r_var = res.res_var
    return m_res, b_res, r_var