#include "m_pd.h"

static t_class *tabredraw_class;

typedef struct _tabredraw {
    t_object x_obj;
    t_symbol *x_arrayname;
} t_tabredraw;

static void tabredraw_bang(t_tabredraw *x)
{
    t_garray *array;
    
    if (!(array = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))) {
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
        return;
    }
    garray_redraw(array);
}

static void *tabredraw_new(t_symbol *s)
{
    t_tabredraw *x = (t_tabredraw *)pd_new(tabredraw_class);
    x->x_arrayname = s;
    return (x);
}

void tabredraw_setup(void)
{
    tabredraw_class = class_new(gensym("tabredraw"),
        (t_newmethod)tabredraw_new,
        0,  // No free needed since we only store a symbol
        sizeof(t_tabredraw),
        0,  // No flags needed
        A_DEFSYM, 0);
        
    class_addbang(tabredraw_class, (t_method)tabredraw_bang);
}
