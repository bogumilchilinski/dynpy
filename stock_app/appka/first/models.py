from django.db import models

# Create your models here.

IDEA_STATUS =(
    ('pending', 'Wating for review'),
    ('accepted', 'Accepted'),
    ('done', 'Done'),
    ('rejected', 'Rejected')
)


class Idea(models.Model):
    title = models.CharField(max_length=255)
    description = models.TextField() 
    status = models.CharField(choices=STATUS, default='pending')


class Vote(models.Model):
    idea = models.ForeignKey(Idea, on_delete = models.CASCADE)
    reason = models.TextField()


