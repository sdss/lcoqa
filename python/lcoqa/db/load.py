#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Oct 15, 2017
# @Filename: load.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import pathlib
import re

import click
import peewee

import playhouse.migrate
import playhouse.reflection

import astropy.io.fits
import astropy.time
import astropy.table

from astropy.utils.console import ProgressBar

from .models import database, Frame, Header, BinTable


type_field_mapping = {int: peewee.IntegerField,
                      float: peewee.FloatField,
                      str: peewee.CharField,
                      bool: peewee.BooleanField}


def get_mjd(date_obs):
    """Returns the LCO SJD."""

    return int(astropy.time.Time(date_obs).mjd + 0.4)


def add_columns(database, table_name, payload={}):
    """Adss new columns to a table, in run-time.

    Parameter:
        database (`Database <http://docs.peewee-orm.com/en/latest/peewee/api.html#Database>`_):
            The peewee database to use.
        table_name (str):
            The name of the table.
        payload (dict):
            A mapping of key-value. All keys that are not present in the table
            will be added. The value will be used to determine the type of the
            new column.

    """

    migrator = playhouse.migrate.SqliteMigrator(database)

    colnames = [col.name for col in database.get_columns(table_name)]

    new_columns = []
    for key in payload:

        db_key = key
        value = payload[key]

        if db_key not in colnames:
            new_field = type_field_mapping[type(value)](null=True)
            new_columns.append(migrator.add_column(table_name, db_key, new_field))

    if len(new_columns) > 0:
        playhouse.migrate.migrate(*new_columns)
        return len(new_columns)
    else:
        return 0


def load_image(image, db, header=True, bintable=True):
    """Loads a proc-gimg to the DB.

    Parameters:
        image (str or `~pathlib.Path`):
            The path to the image to load.
        db (str or `Database <http://docs.peewee-orm.com/en/latest/peewee/api.html#Database>`_):
            A string to the path of a sqlite database (if the path does not
            exists, a new DB will be created) or a pewee database object.
        header (bool):
            Load primary header?
        bintable (bool):
            Load bintable?

    """

    # Loads the image
    image_path = pathlib.Path(image)
    assert image_path.exists()
    assert 'proc-gimg' in image_path.name, 'not a proc-gimg image.'

    header = astropy.io.fits.getheader(image_path, 0)

    if header['IMAGETYP'] == 'dark':
        bintable = None:
    else:
        bintable = astropy.table.Table.read(image_path)

    # Makes sure the DB the models are linked to is pointing to the right DB.

    assert isinstance(db, (str, peewee.SqliteDatabase)), 'invalid db type.'

    if isinstance(db, peewee.SqliteDatabase):
        db = db.database

    if database.database != db:
        if not database.is_closed():
            database.close()
        database.init(db)

    if database.is_closed():
        database.connect()

    # Inserts the tables. If they already exist they won't be replaced.
    database.create_tables([Frame, Header, BinTable], safe=True)

    frame = int(re.search('proc-gimg-([0-9]*)', str(image)).group(1))
    mjd = get_mjd(header['DATE-OBS'])

    # Adds or retrieves the new frame.
    frame_dbo, __ = Frame.get_or_create(frame=frame, mjd=mjd, processed=True)

    # Adds new columns to Header
    header_values = {key.lower().replace('-', '_'): header[key]
                     for key in header if key.lower() != 'comment'}
    add_columns(database, 'header', header_values)

    # Gets the new header model from reflection
    models = playhouse.reflection.Introspector.from_database(database).generate_models()
    full_header_model = models['header']

    with database.atomic():
        header_dbo = full_header_model.get_or_create(frame_pk=frame_dbo.pk, extension=0)[0]
        header_blob = bytes(header.tostring(), 'utf-8')
        full_header_model.update(header_blob=header_blob,
                                 **header_values).where(
                                     full_header_model.pk == header_dbo.pk).execute()

    if bintable is None:
        return

    # Loads the bintable
    db_cols = database.get_columns('bintable')

    # Removes previous bintable row for this frame, just in case.
    with database.atomic():
        BinTable.delete().where(BinTable.frame_pk == frame_dbo.pk).execute()

    # Selects only columns that are in the model.
    bintable_db = bintable[[col.name for col in db_cols if col.name in bintable.colnames]]
    bintable_data = [{bintable_db.colnames[ii]: col for ii, col in enumerate(row)}
                     for row in bintable_db]

    [dd.update({'frame_pk': frame_dbo.pk}) for dd in bintable_data]  # Adds frame_pk
    with database.atomic():
        for data_dict in bintable_data:
            BinTable.create(**data_dict)


# CLI interface. Intended to be added as a subparser for bin/lcoqa.

@click.command()
@click.argument('images', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('dbpath', type=click.Path())
@click.option('-r', '--recursive', is_flag=True,
              help='recursively loads all proc-gimg')
def dbload(images, dbpath, recursive=False):
    """Loads proc- gimages to the DB."""

    paths = []

    for image in images:
        image_path = pathlib.Path(image)
        if image_path.is_dir():
            if recursive:
                paths += image_path.rglob('proc-gimg*')
            else:
                paths += image_path.glob('proc-gimg*')
        else:
            paths.append(image_path)

    with ProgressBar(len(paths)) as bar:
        for path in sorted(paths):
            try:
                load_image(path, dbpath)
            except Exception as ee:
                print('failed loading file {}: {}'.format(str(path), ee))
            bar.update()
