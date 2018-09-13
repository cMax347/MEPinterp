module kubo_ahc
	use constants,		only:	dp
	use statistics,		only:	fermi_dirac
	implicit none


	private
	public			::			kubo_ahc_tens	

contains


	pure function kubo_ahc_tens(en_k, Om_kab, eFermi, T_kelvin, unit_vol) result( o_ahc)
		real(dp),		intent(in)		::	en_k(:), eFermi, T_kelvin, unit_vol
		complex(dp),	intent(in)		::	Om_kab(:,:,:,:)
		real(dp)						::	o_ahc(3,3)
		integer							::	n
		!
		o_ahc	=	0.0_dp
		!
		do n = 1, size(en_k)
			o_ahc	=	o_ahc	-	real(Om_kab(:,:,n,n),dp)		*	fermi_dirac(en_k(n),	eFermi, 	T_kelvin) /	unit_vol
		end do
		!
		!
	end function

end module kubo_ahc