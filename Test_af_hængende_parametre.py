

# In[ ]:


from classy import Class

params = {
    'omega_b': 0.0220156,
    'omega_cdm': 0.12215366,
    '100*theta_s': 1.03940682,
    'ln10^{10}A_s': 3.04884021,
    'n_s': 0.9574973,
    'tau_reio': 0.02692163,
    'r': 0.18079943,
    'm_ncdm_interacting': 0.82079422,
    'log10_G_eff_ur': -2.45657156,
    'deg_ncdm_interacting': 2.19480209,

    'output': 'mPk, tCl, lCl, pCl',
    'l_max_scalars': 2500,
    'modes': 's,t',
    'k_pivot': 0.05,
    'N_ncdm_interacting': 1,
    'ncdm_fluid_approximation': 3,
    'N_ur': 0.0
}
cosmo = Class()
cosmo.set(params)
cosmo.compute()





# In[1]:


from classy_nu import Class

params = {
    'omega_b': 0.02180326,
    'omega_cdm': 0.11229913,
    '100*theta_s': 1.04025133,
    'ln10^{10}A_s': 3.10189877,
    'n_s': 0.97736789,
    'tau_reio': 0.0747002,
    'r': 0.05117754,
    'm_ncdm_interacting': 0.24910741,
    'log10_G_eff_ur': -3.72876086,
    'deg_ncdm_interacting': 3.35010057,

    'output': 'mPk, tCl, lCl, pCl',
    'l_max_scalars': 2500,
    'modes': 's,t',
    'k_pivot': 0.05,
    'N_ncdm_interacting': 1,
    'ncdm_fluid_approximation': 3,
    'N_ur': 0.0
}
cosmo = Class()
cosmo.set(params)
cosmo.compute()



# In[ ]:


from classy_nu import Class

params = {
    'omega_b': 0.02213333,
    'omega_cdm': 0.12016931,
    '100*theta_s': 1.04004526,
    'ln10^{10}A_s': 3.10986326,
    'n_s': 0.97561151,
    'tau_reio': 0.05240832,
    'r': 0.17031911,
    'm_ncdm_interacting': 0.91603865,
    'log10_G_eff_ur': -5.01390696,
    'deg_ncdm_interacting': 2.31943228,
    
    'output': 'mPk, tCl, lCl, pCl',
    'l_max_scalars': 2500,
    'modes': 's,t',
    'k_pivot': 0.05,
    'N_ncdm_interacting': 1,
    'ncdm_fluid_approximation': 3,
    'N_ur': 0.0
}
cosmo = Class()
cosmo.set(params)
cosmo.compute()


# In[ ]:




