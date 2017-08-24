#include "standalone.h"

void setup_async_file(t_standalone *sim)
{
	//makes a file of appropriate size to write the output to.
	char *filename;
	asprintf(&filename, "%s-%d.jgrav", sim->name, sim->ticks_done);
	pthread_mutex_lock(&sim->output_mutex);
	sim->fp = fopen(filename, "w");
	fwrite(&sim->dataset->particle_cnt, sizeof(long), 1, sim->fp);
	fseek(sim->fp, sizeof(cl_float4) * sim->dataset->particle_cnt, SEEK_SET);
	fputc('\0', sim->fp);
	fseek(sim->fp, 0, SEEK_SET);
	pthread_mutex_unlock(&sim->output_mutex);
	free(filename);
}

void async_save(t_standalone *sim, unsigned long offset, t_workunit *wu)
{
	pthread_mutex_lock(&sim->output_mutex);
	fseek(sim->fp, offset * sizeof(cl_float4) + sizeof(long), SEEK_SET);
	fwrite(wu->N, sizeof(cl_float4), wu->n_count, sim->fp);
	pthread_mutex_unlock(&sim->output_mutex);
}