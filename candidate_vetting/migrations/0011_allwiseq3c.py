from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('candidate_vetting', '0010_delvedr3q3c_delete_usergalaxyq3c'),
    ]

    operations = [
        migrations.CreateModel(
            name='AllwiseQ3C',
            fields=[
                ('cntr', models.BigIntegerField(primary_key=True, serialize=False)),
                ('designation', models.CharField(blank=True, max_length=20, null=True)),
                ('ra', models.FloatField(blank=True, null=True)),
                ('dec', models.FloatField(blank=True, null=True)),
                ('w1mpro', models.FloatField(blank=True, null=True)),
                ('w1sigmpro', models.FloatField(blank=True, null=True)),
                ('w1snr', models.FloatField(blank=True, null=True)),
                ('w2mpro', models.FloatField(blank=True, null=True)),
                ('w2sigmpro', models.FloatField(blank=True, null=True)),
                ('w2snr', models.FloatField(blank=True, null=True)),
                ('w3mpro', models.FloatField(blank=True, null=True)),
                ('w3sigmpro', models.FloatField(blank=True, null=True)),
                ('w3snr', models.FloatField(blank=True, null=True)),
                ('w4mpro', models.FloatField(blank=True, null=True)),
                ('w4sigmpro', models.FloatField(blank=True, null=True)),
                ('w4snr', models.FloatField(blank=True, null=True)),
                ('cc_flags', models.CharField(blank=True, max_length=4, null=True)),
                ('ext_flg', models.IntegerField(blank=True, null=True)),
                ('ph_qual', models.CharField(blank=True, max_length=4, null=True)),
            ],
            options={
                'db_table': 'allwise_q3c',
                'managed': False,
            },
        ),
    ]
