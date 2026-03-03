/* beng207_2026-window.c
 *
 * Copyright 2026 shh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "config.h"

#include "beng207_2026-window.h"

struct _Beng2072026Window
{
	AdwApplicationWindow  parent_instance;

	/* Template widgets */
	GtkLabel            *label;
};

G_DEFINE_FINAL_TYPE (Beng2072026Window, beng207_2026_window, ADW_TYPE_APPLICATION_WINDOW)

static void
beng207_2026_window_class_init (Beng2072026WindowClass *klass)
{
	GtkWidgetClass *widget_class = GTK_WIDGET_CLASS (klass);

	gtk_widget_class_set_template_from_resource (widget_class, "/org/gnome/Example/beng207_2026-window.ui");
	gtk_widget_class_bind_template_child (widget_class, Beng2072026Window, label);
}

static void
beng207_2026_window_init (Beng2072026Window *self)
{
	gtk_widget_init_template (GTK_WIDGET (self));
}
