#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <gdk/gdkkeysyms.h>
#include "genometools.h"
#include "core/bioseq.h"

GtDiagram *d = NULL;
GtLayout *l = NULL;
GtStyle *sty = NULL;
GtError *err = NULL;
GtkWidget *area;

static gboolean on_expose_event(GtkWidget *widget,
                                GdkEventExpose *event,
                                gpointer data)
{
  cairo_t *cr;
  GtCanvas *canvas = NULL;
  unsigned long height;
  int rval;
  if (!d || widget->allocation.width <= 30) return FALSE;

  /* render image */
  l = gt_layout_new(d, widget->allocation.width, sty, err);
  if (!l) return FALSE;
  rval = gt_layout_get_height(l, &height, err);
  gt_assert(rval == 0);
  gtk_layout_set_size(GTK_LAYOUT(widget),
                      widget->allocation.width,
                      height);
  cr = gdk_cairo_create(GTK_LAYOUT(widget)->bin_window);
  cairo_rectangle(cr, event->area.x, event->area.y, event->area.width,
                  event->area.height);
  cairo_clip(cr);
  canvas = gt_canvas_cairo_context_new(sty, cr, 0, widget->allocation.width,
                                       height, NULL, err);
  gt_assert(canvas);
  gt_layout_sketch(l, canvas, err);
  gt_layout_delete(l);
  gt_canvas_delete(canvas);
  return FALSE;
}

static void
open_file(GtkWidget *widget,  gpointer user_data)
{
  GtkWidget *dialog;
  GtkFileFilter *gff3filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(gff3filter, "*.gff3");
  dialog = gtk_file_chooser_dialog_new ("Open File",
                GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(widget))),
                GTK_FILE_CHOOSER_ACTION_OPEN,
                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                NULL);
  gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), gff3filter);

  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    const char *seqid;
    GtGenomeNode *gn = NULL;
    GtFeatureIndex *features = NULL;
    GtRange qry_range;
    int had_err = 0;
    gt_error_unset(err);
    /* file given, load GFF file */
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    features = gt_feature_index_memory_new();
    gt_feature_index_add_gff3file(features,
                                  filename,
                                  err);

    if (!gt_error_is_set(err))
    {
      GtkWidget *w = GTK_WIDGET(area);
      GtkListStore *store;
      GtBioseq *bioseq;
      seqid = gt_feature_index_get_first_seqid(features, err);
      gt_feature_index_get_range_for_seqid(features, &qry_range, seqid, err);
      gt_diagram_delete(d);
      d = gt_diagram_new(features, seqid, &qry_range, sty, err);
      gtk_widget_queue_draw_area(w, 0, 0, w->allocation.width,
                                 w->allocation.height);
      gtk_widget_destroy(dialog);
    } else {
      GtkWidget *edialog;
      gtk_widget_destroy(dialog);
      edialog =  gtk_message_dialog_new(GTK_WINDOW(
                                          gtk_widget_get_toplevel(
                                             GTK_WIDGET(widget))),
                                        0,
                                        GTK_MESSAGE_ERROR,
                                        GTK_BUTTONS_OK,
                                        "Error loading file '%s' : %s",
                                        filename, gt_error_get(err));
      gtk_dialog_run(GTK_DIALOG(edialog));
      gtk_widget_destroy (edialog);
    }
    g_free(filename);
  } else {
    gtk_widget_destroy(dialog);
  }
}

int main(int argc, char *argv[])
{
  GtkWidget *window, *vbox, *sw;
  GtkWidget *menubar = NULL;
  GtkWidget *filemenu = NULL;
  GtkWidget *open = NULL, *openbt = NULL;
  GtkWidget *file = NULL;
  GtkWidget *quit = NULL;
  GtkWidget *sep = NULL;
  GtkWidget *vp = NULL;
  GtkWidget *sbar = NULL;
  GtkWidget *tb = NULL;
  GtkAdjustment *vadj = NULL;
  GtkAccelGroup *accel_group = NULL;
  GtkCellRenderer *cell;

  gtk_set_locale();

  gtk_init(&argc, &argv);

  gt_lib_init();
  err = gt_error_new();
  sty = gt_style_new(err);
  gt_style_load_file(sty, "gtdata/sketch/default.style", err);

  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_default_size (GTK_WINDOW (window), 700, 598);
  gtk_window_set_title(GTK_WINDOW(window), "AnnotationSketch viewer demo");
  vbox = gtk_vbox_new(FALSE, 0);
  sw = GTK_WIDGET(gtk_scrolled_window_new(NULL, vadj));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),GTK_POLICY_NEVER,
                                                         GTK_POLICY_AUTOMATIC);
  tb = gtk_toolbar_new();
  openbt =  GTK_WIDGET(gtk_tool_button_new_from_stock(GTK_STOCK_OPEN));
  gtk_toolbar_insert(GTK_TOOLBAR(tb), GTK_TOOL_ITEM(openbt), -1);

  cell = gtk_cell_renderer_text_new();
  sbar = gtk_statusbar_new();

  area = gtk_layout_new(NULL, vadj);
  gtk_layout_set_size(GTK_LAYOUT(area), 600, 500);
  gtk_container_add(GTK_CONTAINER(sw), area);

  gtk_container_add(GTK_CONTAINER(window), vbox);
  menubar = gtk_menu_bar_new();
  filemenu = gtk_menu_new();

  accel_group = gtk_accel_group_new();
  gtk_window_add_accel_group(GTK_WINDOW(window), accel_group);

  file = gtk_menu_item_new_with_mnemonic("_File");
  open = gtk_image_menu_item_new_from_stock(GTK_STOCK_OPEN, NULL);
  sep = gtk_separator_menu_item_new();
  quit = gtk_image_menu_item_new_from_stock(GTK_STOCK_QUIT, accel_group);

  gtk_widget_add_accelerator(quit, "activate", accel_group,
      GDK_q, GDK_CONTROL_MASK, GTK_ACCEL_VISIBLE);
  gtk_widget_add_accelerator(open, "activate", accel_group,
      GDK_o, GDK_CONTROL_MASK, GTK_ACCEL_VISIBLE);

  gtk_menu_item_set_submenu(GTK_MENU_ITEM(file), filemenu);
  gtk_menu_shell_append(GTK_MENU_SHELL(filemenu), open);
  gtk_menu_shell_append(GTK_MENU_SHELL(filemenu), sep);
  gtk_menu_shell_append(GTK_MENU_SHELL(filemenu), quit);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), file);
  gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), tb, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), sbar, FALSE, FALSE, 0);

  g_signal_connect(area, "expose-event",
      G_CALLBACK (on_expose_event), NULL);
  g_signal_connect(window, "destroy",
      G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect(openbt, "clicked",
      G_CALLBACK (open_file), NULL);
  g_signal_connect(quit, "activate",
      G_CALLBACK (gtk_main_quit), NULL);

  gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
  gtk_widget_show_all(window);
  gtk_main();

  gt_error_delete(err);
  gt_style_delete(sty);
  gt_diagram_delete(d);
  return 0;
}
