/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   calc_thread.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: cyildiri <cyildiri@student.42.us.org>      +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/02 18:29:42 by cyildiri          #+#    #+#             */
/*   Updated: 2017/06/15 01:54:57 by cyildiri         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "standalone.h"

typedef struct s_calc
{
	t_standalone * sim;
	int cq_id;
}				t_calc;

static void *calc_thread(void *param)
{
	t_standalone *sim;
	t_workunit  *workunit;

	sim = ((t_calc *)param)->sim;
	int cq_id = ((t_calc *)param)->cq_id;
	while (1)
	{
		if (sem_wait(sim->calc_thread_sem) < 0)
		    printf("CALC- sem_wait failed with err:%d\n", errno);
		workunit = queue_pop(&sim->todo_work);
		if(!workunit)
			printf("WU was NULL!!\n");
		do_workunit(sim, workunit, cq_id);
		queue_enqueue(&sim->completed_work, queue_create_new(workunit));
    	sem_post(sim->integration_thread_sem);
	}
	return (0);
}

void launch_calculation_thread(t_standalone *sim)
{
	//this should launch GPU_COUNT threads
	pthread_t *threads = (pthread_t *)calloc(sim->context->device_count, sizeof(pthread_t));
	for (int i = 0; i < sim->context->device_count; i++)
	{
		t_calc *calc = calloc(1, sizeof(t_calc));
		calc->sim = sim;
		calc->cq_id = i;
		pthread_create(&(threads[i]), NULL, calc_thread, calc);
	}
}
