								for (sp2=(sp+1); sp2<=standingdiv; ++sp2)	{									species2=extants[sp2];									if (species2>observed[sampled-1])	sp2=standingdiv;									else	{										for (f=0; f<sampled; ++f)	{											if (observed[f]==species2 && tex[species2]>stage)	{												mod=dmin((tex[species2]-stage),1.0);												v=1-pow(1-(fr*vfr),mod);												find=(((unsigned int) rand())%RAND_MAX)/((double) RAND_MAX);												if (find <= v)	la[f]=stage;												f=sampled;												}	/* a previously sampled & extant species is found	*/											}	/* go through sampled species	*/										}	/* if this might be a found species, check to see if it is	*/									}		/* */