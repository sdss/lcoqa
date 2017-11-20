#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Nov 2, 2017
# @Filename: scaleml.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import pathlib

import astropy.table as table

import playhouse.reflection
from peewee import fn

from ..db.models import database


def get_data(guider_db_file, min_guiderms=0.4):
    """Returns a table with the best frame per plate observed."""

    results_data = []

    guider_db_path = pathlib.Path(guider_db_file)
    assert guider_db_path.exists()

    database.init(str(guider_db_path))
    models = playhouse.reflection.Introspector.from_database(database).generate_models()

    Frame = models['frame']
    Header = models['header']

    plateid_unique = Header.select(Header.plateid).distinct().where(
        Header.plateid.is_null(False)).tuples()

    for plateid in plateid_unique:

        headers = Header.select().where(Header.plateid == plateid[0])
        mjds = headers.join(Frame).select(Frame.mjd).distinct().tuples()

        for mjd in mjds:
            headers_mjd = headers.join(Frame).where(Frame.mjd == mjd)
            valid_headers = headers_mjd.where(Header.imagetyp == 'object').where(
                Header.trpos != 'NaN').where((Header.guideaxe == 1) & (Header.guidefoc == 1))

            min_guiderms = valid_headers.select(fn.MIN(Header.gdrms)).where(Header.gdrms > 0)
            min_guiderms_headers = valid_headers.where(Header.gdrms == min_guiderms)

            if min_guiderms_headers.count() != 1:
                continue

            row = min_guiderms_headers[0]

            results_data.append((row.plateid, row.frame_pk.mjd, row.ra, row.dec, row.t_out,
                                 row.t_in, row.t_prim, row.t_cell, row.t_floor, row.trpos,
                                 row.gdrms, row.armass, row.m2piston))

    results_table = table.Table(rows=results_data,
                                names=['plateifu', 'mjd', 'ra', 'dec', 't_out', 't_in',
                                       't_prim', 't_cell', 't_floor', 'trpos', 'gdrms',
                                       'armass', 'm2piston'])

    return results_table
