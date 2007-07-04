/*
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DIAGRAM_H
#define DIAGRAM_H

#include <libgtcore/array.h>
#include <libgtview/config.h>
#include <libgtcore/range.h>
#include <libgtcore/hashtable.h>

typedef struct Diagram Diagram;

             /*!
            Initialize a new diagram object.
            \param features pointer to the array of genome nodes.
            \param range the given range of the diagam.
            \param config pointer to the configuration object.
            \param env Pointer to Environment object.
            */
Diagram*    diagram_new(Array* features,Range,Config*,Env*);

            /*!
            delivers the range of a diagram
            \param diagram Pointer to diagram object.
            */
Range       diagram_get_range(Diagram* diagram);

            /*!
            Update the configuration object with new settings.
            \param diagram pointer to the diagram object.
            \param config pointer to the configuration object.
            \param env Pointer to Environment object.
            */
void        diagram_set_config(Diagram*,Config*,Env*);

            /*!
            Delivers the hashtable with the stored tracks.
            \param diagram pointer to the diagram object.
            */
Hashtable*  diagram_get_tracks(Diagram*);

            /*!
            Returns the number of all lines in the diagram.
            \param diagram pointer to the diagram object.
            \param env Pointer to Environment object.
            */
int         diagram_get_total_lines(Diagram*, Env*);

            /*!
            Returns the number of all lines in the diagram.
            \param diagram pointer to the diagram object.
            \param env Pointer to Environment object.
            */
int         diagram_get_number_of_tracks(Diagram *diagram);

            /*
            Delete the diagram object.
            \param diagram pointer to the diagram object.
            \param env Pointer to Environment object.
            */
void        diagram_delete(Diagram*,Env*);

            /*
            generate a feature index test structure and
            test the diagram functions.
            \param env Pointer to Environment object.
            */
int         diagram_unit_test(Env*);

#endif
